import logging
import sys

import torch
import pandas as pd
import os

import vcuppa_train

logger = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s %(levelname)5s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

# During training, we use a specific Key ordering. We must use the same.
def load_cuppa_features(cuppa_features_path: str, definitions: pd.DataFrame):
    cuppa_feature_df = pd.read_csv(cuppa_features_path, sep="\t")
    key_order = definitions["Key"].to_list()
    key_index = cuppa_feature_df.columns.to_list().index("Key")
    cuppa_feature_df = cuppa_feature_df.set_index("Key").reindex(key_order).reset_index()
    cuppa_feature_df.insert(key_index, "Key", cuppa_feature_df.pop("Key"))
    return cuppa_feature_df.reset_index(drop=True)

def predict_sample1(df, driver_df, snv_df, purple_plots_dir, model):

    df = df[[c for c in df.columns if c in ["sampleId", "wgsSampleId", "cancerSubtype"]]]
    df = vcuppa_train.filter_df(df).reset_index(drop=True)
    snv_df = snv_df[snv_df["Category"] == "sig"]

    # remove any sample that is not driver df
    df = df[df["sampleId"].isin(driver_df.columns)].reset_index(drop=True)

    # add all the driver to the df
    df["drivers"] = df["sampleId"].apply(lambda x: driver_df[x].to_numpy())
    df["snvs"] = df["sampleId"].apply(lambda x: snv_df[x].to_numpy())

    model.eval()
    cnn_pred_score = []
    cnn_pred = []

    i = 0

    with torch.no_grad():
        for idx, row in df.iterrows():
            sample_id = row["sampleId"]
            circos_png_path = f'{purple_plots_dir}/{sample_id}/{sample_id}.circos.png'

            if not os.path.exists(circos_png_path):
                print(f"{sample_id} has no purple plot")

                cnn_pred_score.append(float('nan'))
                cnn_pred.append(float('nan'))

                continue

            image_tensor = vcuppa_train.image_to_tensor(circos_png_path).unsqueeze(0)
            driver_tensor = vcuppa_train.driver_to_tensor(row["drivers"], row["snvs"]).unsqueeze(0)
            pred_tensor = model(image_tensor, driver_tensor).squeeze(0)
            pred_tensor = torch.nn.functional.softmax(pred_tensor, dim=0)
            pred_score = pred_tensor[torch.argmax(pred_tensor, dim=0)].item()
            pred = vcuppa_train.tensor_to_cancer_type(pred_tensor)

            # print(f"pred: {pred}, tensor: {pred_tensor}")

            cnn_pred_score.append(pred_score)
            cnn_pred.append(pred)
            i += 1

            if (i % 100) == 0:
                print(f"{i} complete")

    df["cnnPred"] = cnn_pred
    df["cnnPredScore"] = cnn_pred_score

    df["cnnCorrect"] = (df["cancerSubtype"] == df["cnnPred"])

    return df


def predict_sample(sample, purple_circos_path, cuppa_features, model):
    # Load necessary dataframes and directories

    image_tensor = vcuppa_train.image_to_tensor(purple_circos_path).unsqueeze(0)

    # get the drivers and sig out of the cuppa_features
    drivers = cuppa_features[cuppa_features["Category"] != "sig"]["Value"]
    sigs = cuppa_features[cuppa_features["Category"] == "sig"]["Value"]

    driver_tensor = vcuppa_train.driver_to_tensor(drivers, sigs).unsqueeze(0)
    pred_tensor = model(image_tensor, driver_tensor).squeeze(0)
    pred_tensor = torch.nn.functional.softmax(pred_tensor, dim=0)
    pred_score = pred_tensor[torch.argmax(pred_tensor, dim=0)].item()
    pred = vcuppa_train.tensor_to_cancer_type(pred_tensor)

    logger.info(f"pred: {pred}, tensor: {pred_tensor}, pred_score: {pred_score}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="run vcuppa prediction on a sample")
    parser.add_argument('--sample', help='sample Id', required=True)
    parser.add_argument('--purple_circos', help='path to purple circos', required=True)
    parser.add_argument('--cuppa_data', help='path to cuppa sample data', required=True)
    parser.add_argument('--model', help='path to torch script model', required=True)
    parser.add_argument('--cuppa_data_definitions', help='path to cuppa data definitions', required=True)
    parser.add_argument('--output_dir', help='path to output directory', required=True)
    args = parser.parse_args()

    model = torch.jit.load(args.model, map_location=torch.device('cpu'))
    model.eval()
    cuppa_features_def_df = pd.read_csv(args.cuppa_feature_definitions, sep="\t")
    cuppa_features_df = load_cuppa_features(args.cuppa_features, cuppa_features_def_df)

    with torch.no_grad():
        predict_sample(args.sample, args.purple_circos, cuppa_features_df, model)


if __name__ == "__main__":
    main()
