"""Pathology foundation models.

Wraps the UNI2-h histopathology foundation model (loaded through ``timm``),
applies the model-specific preprocessing (resize to 224x224 and ImageNet
normalisation) and extracts patch embeddings in batches.
"""

import os
from typing import Optional

import numpy as np
import timm
import torch
import torchvision.transforms as transforms
from PIL import Image
from dotenv import load_dotenv
from tqdm import tqdm

# Load environment variables (such as HF_TOKEN) from a .env file if present
load_dotenv()

UNI2_MODEL_NAME = "hf_hub:MahmoodLab/UNI2-h"
UNI2_INPUT_SIZE = 224
UNI2_EMBED_DIM = 1536
IMAGENET_MEAN = (0.485, 0.456, 0.406)
IMAGENET_STD = (0.229, 0.224, 0.225)


def _tqdm_position() -> int:
    """Progress-bar row, so parallel workers do not overwrite each other's bars."""
    return int(os.environ.get("TQDM_POSITION", "0"))


def get_uni2_kwargs() -> dict:
    """Architecture keyword arguments required to instantiate UNI2-h.

    These reproduce the ViT-Giant/SwiGLU configuration published by the UNI2
    authors; changing any of them yields a model whose weights will not load.
    """
    return {
        "img_size": UNI2_INPUT_SIZE,
        "patch_size": 14,
        "depth": 24,
        "num_heads": 24,
        "init_values": 1e-5,
        "embed_dim": UNI2_EMBED_DIM,
        "mlp_ratio": 2.66667 * 2,
        "num_classes": 0,                       # feature extractor: no classification head
        "no_embed_class": True,
        "mlp_layer": timm.layers.SwiGLUPacked,
        "act_layer": torch.nn.SiLU,             # SiLU is standard for SwiGLU architectures
        "reg_tokens": 8,
        "dynamic_img_size": True,
    }


class UNI2FeatureExtractor:
    """Feature extractor for the UNI2-h pathology foundation model.

    Loads the ViT-Giant/SwiGLU model (``embed_dim`` = 1536) and runs batched
    inference over patch arrays.

    Weights come from the Hugging Face Hub by default, which requires accepting
    the model licence and being authenticated (``HF_TOKEN``); pass
    ``weights_path`` to load a local checkpoint instead.
    """

    def __init__(self, weights_path: Optional[str] = None, device: Optional[str] = None):
        """
        Args:
            weights_path: Path to a local PyTorch checkpoint (``.pth`` / ``.pt``).
                If ``None``, the weights are downloaded from the Hugging Face Hub.
            device: Target device (``'cuda'``, ``'cpu'``). Auto-detected if ``None``.

        Raises:
            FileNotFoundError: If ``weights_path`` is given but does not exist.
        """
        self.device = device if device else ("cuda" if torch.cuda.is_available() else "cpu")
        print(f"Initializing UNI2 (UNI2-h) on device: {self.device}")

        timm_kwargs = get_uni2_kwargs()
        pretrained = weights_path is None
        self.model = timm.create_model(UNI2_MODEL_NAME, pretrained=pretrained, **timm_kwargs)

        if weights_path:
            if not os.path.exists(weights_path):
                raise FileNotFoundError(f"Local weights path not found: {weights_path}")
            print(f"Loading custom UNI2 weights from checkpoint: {weights_path}")
            state_dict = torch.load(weights_path, map_location="cpu")
            if "state_dict" in state_dict:
                state_dict = state_dict["state_dict"]
            self.model.load_state_dict(state_dict)

        self.model = self.model.to(self.device)
        self.model.eval()

        # Standard UNI2 preprocessing: resize the tile to 224x224, then normalise
        # with the ImageNet statistics the model was trained with.
        self.transform = transforms.Compose([
            transforms.Resize((UNI2_INPUT_SIZE, UNI2_INPUT_SIZE)),
            transforms.ToTensor(),
            transforms.Normalize(mean=IMAGENET_MEAN, std=IMAGENET_STD),
        ])

    @torch.no_grad()
    def extract_features(self, patches: np.ndarray, batch_size: int = 64) -> np.ndarray:
        """Extract embeddings from a collection of raw patches.

        Args:
            patches: Array of shape (P, H, W, 3) of uint8.
            batch_size: Patches per forward pass; lower it if the GPU runs out of memory.

        Returns:
            Embedding matrix of shape (P, embed_dim) of float32. An empty input
            yields an empty array with the correct embedding dimension.
        """
        num_patches = patches.shape[0]
        if num_patches == 0:
            # Probe the model once to report the right embedding dimension
            dummy_input = torch.zeros((1, 3, UNI2_INPUT_SIZE, UNI2_INPUT_SIZE), device=self.device)
            embed_dim = self.model(dummy_input).shape[1]
            return np.empty((0, embed_dim), dtype=np.float32)

        embeddings_list = []
        for idx in tqdm(range(0, num_patches, batch_size), desc="UNI2 Inference",
                        unit="batch", leave=False, position=_tqdm_position()):
            batch_patches_np = patches[idx: idx + batch_size]

            transformed = [self.transform(Image.fromarray(patch)) for patch in batch_patches_np]
            batch_tensor = torch.stack(transformed).to(self.device)

            batch_embeddings = self.model(batch_tensor)
            embeddings_list.append(batch_embeddings.cpu().numpy())

        return np.concatenate(embeddings_list, axis=0).astype(np.float32)
