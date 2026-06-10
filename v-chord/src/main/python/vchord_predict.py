import os
import pandas as pd
import torch
import torch.utils.data as data_utils
from torchvision.io import read_image
import vchord_train
from common import LOGGER, DEVICE


class PredictionDataset(data_utils.Dataset):
    """Dataset for batch prediction"""

    def __init__(self, df: pd.DataFrame, purple_root: str, transform) -> None:
        self.df = df
        self.purple_root = purple_root
        self.transform = transform

    def __len__(self) -> int:
        return len(self.df)

    def __getitem__(self, idx: int) -> tuple[torch.Tensor, torch.Tensor, int]:
        row = self.df.iloc[idx]
        sample_id = row["sampleId"]
        image_path = f'{self.purple_root}/{sample_id}.circos.png'

        # Check if image exists
        if not os.path.exists(image_path):
            # Return empty tensors if image doesn't exist
            image_tensor = torch.zeros(3, vchord_train.IMAGE_SIZE, vchord_train.IMAGE_SIZE)
            type_tensor = torch.zeros(vchord_train.NUM_CANCER_TYPES)
            return image_tensor, type_tensor, -1  # -1 indicates missing image

        # Load and transform image
        image = read_image(image_path)
        image_tensor = self.transform(image)

        # Get cancer type and purity
        cancer_type = row.get("primaryTumorLocation", None)
        purity = row.get("purity", 1.0)
        type_tensor = vchord_train.cancer_type_to_tensor(cancer_type, purity)

        return image_tensor, type_tensor, idx


def predict_batch(model: torch.nn.Module, dataloader: data_utils.DataLoader, num_samples: int) -> list[float]:
    """Run batch predictions and return list of predictions"""
    model.eval()

    predictions = [float('nan')] * num_samples
    num_processed = 0
    num_batches = len(dataloader)

    with torch.no_grad():
        for batch_idx, (image_batch, type_batch, idx_batch) in enumerate(dataloader):
            # Move to device
            image_batch = image_batch.to(DEVICE)
            type_batch = type_batch.to(DEVICE)

            # Get predictions and apply sigmoid
            logits = model(image_batch, type_batch)
            probs = torch.sigmoid(logits).cpu().numpy().flatten()

            # Store predictions at correct indices
            for i, idx in enumerate(idx_batch):
                if idx != -1:  # Only store if image existed
                    predictions[idx] = float(probs[i])

            num_processed += len(idx_batch)
            LOGGER.info(f"Processed batch {batch_idx + 1}/{num_batches}, {num_processed}/{num_samples} samples")

    return predictions


def load_model(model_path: str) -> torch.nn.Module:
    """Load trained model from file"""
    LOGGER.info(f"Loading model from {model_path}")
    model = torch.jit.load(model_path, map_location=torch.device(DEVICE))
    model = model.to(DEVICE)
    model.eval()
    return model


def add_predictions_to_df(df: pd.DataFrame, predictions: list[float]) -> pd.DataFrame:
    """Add prediction columns to dataframe"""
    df["cnnPred"] = [round(p, 3) if not pd.isna(p) else float('nan') for p in predictions]

    # Classify as HR deficient or proficient
    df["cnnHrStatus"] = "HR_PROFICIENT"
    df.loc[df["cnnPred"] > 0.5, "cnnHrStatus"] = "HR_DEFICIENT"
    df.loc[df["cnnPred"].isna(), "cnnHrStatus"] = None

    return df


def log_prediction_summary(df: pd.DataFrame) -> None:
    """Log summary statistics of predictions"""
    num_predicted = df["cnnPred"].notna().sum()
    num_missing = df["cnnPred"].isna().sum()
    num_hrd = (df["cnnPred"] > 0.5).sum()
    num_hrp = (df["cnnPred"] <= 0.5).sum()

    LOGGER.info(f"Prediction summary:")
    LOGGER.info(f"  Predicted: {num_predicted}/{len(df)}")
    LOGGER.info(f"  Missing images: {num_missing}")
    LOGGER.info(f"  HR_DEFICIENT: {num_hrd}")
    LOGGER.info(f"  HR_PROFICIENT: {num_hrp}")


def predict_main(sample_tsv: str, purple_root: str, model_path: str, output_tsv: str,
                 batch_size: int, train_set_tsv: str = None) -> None:
    """
    Run predictions on samples using trained model

    Args:
        sample_tsv: Path to TSV file with sample IDs and metadata
        purple_root: Root directory containing circos PNG files
        model_path: Path to saved model (.pth file)
        output_tsv: Path for output TSV with predictions
        batch_size: Batch size for prediction
        train_set_tsv: Optional path to train_set.tsv.gz to mark training samples
    """
    # Load samples
    LOGGER.info(f"Loading samples from {sample_tsv}")
    df = pd.read_csv(sample_tsv, sep='\t')
    LOGGER.info(f"Loaded {len(df)} samples")

    # Load model
    model = load_model(model_path)

    # Create dataset and dataloader
    transform = vchord_train.make_transform()
    dataset = PredictionDataset(df, purple_root, transform)
    dataloader = data_utils.DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)

    # Run predictions
    LOGGER.info(f"Running predictions with batch_size={batch_size} on {DEVICE}")
    predictions = predict_batch(model, dataloader, len(df))

    # Process predictions
    df = add_predictions_to_df(df, predictions)

    # Merge training set info if provided
    if train_set_tsv:
        LOGGER.info(f"Loading training set from {train_set_tsv}")
        train_df = pd.read_csv(train_set_tsv, sep='\t')
        df["inTrainingSet"] = df["sampleId"].isin(train_df["sampleId"])

    # Save results
    LOGGER.info(f"Saving predictions to {output_tsv}")
    df.to_csv(output_tsv, sep='\t', index=False)

    # Log summary
    log_prediction_summary(df)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Run HRD predictions on samples")
    parser.add_argument('--sample_tsv', help='Input TSV file with sample IDs', required=True)
    parser.add_argument('--purple_root', help='Root directory containing circos PNG files', required=True)
    parser.add_argument('--model_path', help='Path to trained model (.pth or .pt file)', required=True)
    parser.add_argument('--output_tsv', help='Output TSV file with predictions', required=True)
    parser.add_argument('--batch_size', help='Batch size for prediction', type=int, default=32)
    parser.add_argument('--train_set_tsv', help='Optional path to train_set.tsv.gz to mark training samples', default=None)
    args = parser.parse_args()

    LOGGER.info(f"Using {DEVICE} device")
    predict_main(args.sample_tsv, args.purple_root, args.model_path, args.output_tsv,
                 args.batch_size, args.train_set_tsv)


if __name__ == "__main__":
    main()
