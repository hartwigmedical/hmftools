import pandas as pd

import numpy as np
import math
import os, sys, time

import torch
from torch import nn
import torch.utils.data as data_utils

import torchvision
from torch.utils.data import DataLoader
# use v2 as it claims to be faster
from torchvision.transforms import v2
from torchvision.io import read_image

import matplotlib.pyplot as plt

import logging

import vcuppa_model

from vcuppa_model import cancer_types
from vcuppa_model import IMAGE_SIZE

logger = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s %(levelname)5s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logger.setLevel(logging.DEBUG)

# -1.0 = 69%, -0.6=72% -0.5 = 72%, -0.4=73% -0.3 = 73%, -0.2=76%, -0.1 = 75%
DRIVER_FILL_NAN = -0.2

# Get cpu or gpu device for training.
device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"


# convert a df to tensor to be used in pytorch
def df_to_tensor(df):
    return torch.from_numpy(df.values).float().to(device)


def filter_df(df):
    if 'qcStatus' in df.columns:
        df = df[df["qcStatus"] == "PASS"]
    return df[df["cancerSubtype"].isin(cancer_types)]


# create a transform object to convert image to pytorch tensor
def make_transform(image_size=IMAGE_SIZE):
    x = [v2.Resize(image_size),
         v2.ToDtype(torch.float32, scale=True)]
    return v2.Compose(x)


def image_to_tensor(circos_png_path, image_size=IMAGE_SIZE, transform=None):
    if transform is None:
        transform = make_transform(image_size)
    circos_png = transform(read_image(circos_png_path))
    return circos_png
    #return torch.empty((IMAGE_SIZE, IMAGE_SIZE))


def cancer_type_to_tensor(cancer_type: str):
    a = [0.0] * len(cancer_types)
    for i in range(len(cancer_types)):
        if cancer_type == cancer_types[i]:
            a[i] = 1.0
    # logger.info(f"cancer type: {cancer_type}, encoded: {a}")
    return torch.tensor(a, dtype=torch.float32)

def tensor_to_cancer_type(tensor: torch.Tensor) -> str:
    assert tensor.size(dim=0) == len(cancer_types), \
        f"Tensor length {tensor.size(dim=0)} does not match the number of cancer types {len(cancer_types)}"
    return cancer_types[torch.argmax(tensor).item()]

def driver_to_tensor(driver: np.ndarray, snv: np.ndarray):
    a = np.concatenate([driver, np.log1p(snv)])
    return torch.from_numpy(a).nan_to_num(nan=DRIVER_FILL_NAN).float()


class VCuppaDataset(data_utils.Dataset):

    def __init__(self, df, augment):
        start = time.time()

        # duplicate all certain HRD samples to enrich the size
        # df = pd.concat([df] + [df[df["hrd"] >= 0.8]] * hrd_sample_dup, ignore_index=True)

        self.targets = []
        self.image_tensors = []
        self.driver_tensors = []

        transform = make_transform(IMAGE_SIZE)

        logger.info(f"start loading images, augment={augment}")
        # load all the images into an array
        for idx, row in df.iterrows():
            if not os.path.exists(row["circosPngPath"]):
                logger.warning(f"sample circos: {row['circosPngPath']} not found, skipping")
                continue
            self.image_tensors.append(image_to_tensor(row["circosPngPath"], IMAGE_SIZE, transform))
            self.driver_tensors.append(driver_to_tensor(row["drivers"], row["snvs"]))
            self.targets.append(cancer_type_to_tensor(row["cancerSubtype"]))

        self.transform = None

        elapsed_sec = int(time.time() - start)
        minute, second = divmod(elapsed_sec, 60)
        logger.info(f"loading dataset of size {len(self)}, took {minute:.0f}m {second:.0f}s")

    def __len__(self):
        return len(self.image_tensors)

    def __getitem__(self, idx):
        # we do the rotation here to support augmentation of data
        image_tensor = self.image_tensors[idx]
        if self.transform:
            # NOTE: this is the right place to put it onto GPU, such that the rotation is done on GPU,
            # otherwise the CPU would struggle with the load and become the bottleneck
            # However if we are not doing transform then we keep it in GPU and move whole batch to GPU
            # in the train / test functions
            image_tensor = image_tensor.to(device)
            image_tensor = self.transform(image_tensor)
        driver_tensor = self.driver_tensors[idx]
        target_tensor = self.targets[idx]
        return image_tensor, driver_tensor, target_tensor

    def get_weights(self):
        # Count occurrences of each cancer type
        cancer_counts = [0] * len(cancer_types)
        for target in self.targets:
            cancer_index = torch.argmax(target).item()
            cancer_counts[cancer_index] += 1

        # Compute weights inversely proportional to count
        nominator = len(self.targets) / len(cancer_types)
        weights = torch.tensor([nominator / count if count > 0 else 1.0 for count in cancer_counts],
                               dtype=torch.float32)

        logger.info(f"cancer type weights = {weights}")
        return weights


class EpochStats:
    def __init__(self):
        self.loss_sum = 0.0
        self.loss_count = 0
        self.count = 0
        self.num_correct = 0

    @property
    def loss(self):
        return self.loss_sum / self.loss_count

    @property
    def accuracy(self):
        return self.num_correct / self.count

    @property
    def correct(self):
        return self.num_correct

    def update_loss(self, loss_val, count):
        self.loss_sum += loss_val * count
        self.loss_count += count

    def update_accuracies(self, pred: torch.Tensor, target: torch.Tensor):
        # logger.debug(f"orig pred: \n{pred}")
        # see https://stackoverflow.com/questions/57208913/is-there-a-way-to-do-pytorch-element-wise-equality-treating-each-dimension-as-an

        # force prediction to 0 and 1, see https://discuss.pytorch.org/t/set-max-value-to-1-others-to-0/44350/7
        idx = torch.argmax(pred, dim=1, keepdims=True)
        pred = torch.zeros_like(pred).scatter(1, idx, 1.)
        self.count += pred.size(dim=0)
        self.num_correct += torch.all(torch.eq(pred, target), dim=1).sum().item()
        # logger.debug(f"calc_error: pred: \n{pred}, target: \n{target}, num correct: {self.num_correct}")

    def log(self, name, epoch):
        logger.info(f"[{name} {epoch:>4d}]    loss: {self.loss:>8.5f}    "
                    f"acc:{self.accuracy * 100:>6.2f}% [{self.correct:>4d}/{self.count:>4d}]    ")
        # f"TP:{self.true_pos_rate * 100:>6.2f}% [{self.true_pos:>4d}/{self.num_pos:>4d}]    "
        # f"TN:{self.true_neg_rate * 100:>6.2f}% [{self.true_neg:>4d}/{self.num_neg:>4d}]")


# from the input pandas dataframe create the train / test dataloaders requied for
# pytorch training.
# return (train_dataloader, test_dataloader)
def create_dataloader(df, batch_size, augment, test_fraction):
    # Using Skicit-learn to split data into training and testing sets
    from sklearn.model_selection import train_test_split

    # Split the data into training and testing sets
    train_df, test_df = train_test_split(df, test_size=test_fraction, random_state=0)

    # write out the train and test set
    train_df[["sampleId"]].to_csv("train_set.tsv.gz", sep="\t", index=False)
    test_df[["sampleId"]].to_csv("test_set.tsv.gz", sep="\t", index=False)

    logger.info(f"train size: {len(train_df)}, test size: {len(test_df)}, batch size: {batch_size}, augment: {augment}")

    # Create data loaders.
    train_dataset = VCuppaDataset(train_df, augment=augment)
    train_dataloader = data_utils.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

    # we do not augment testing data
    test_dataset = VCuppaDataset(test_df, augment=False)
    test_dataloader = data_utils.DataLoader(test_dataset, batch_size=batch_size, shuffle=True)

    # for seq_x, linear_x, y, bq, bqr in test_dataloader:
    #    print(f"Shape of seq x [N, C, H]: {seq_x.shape}")
    #    print(f"Shape of linear x [N, C]: {linear_x.shape}")
    #    print(f"Shape of y: {y.shape} {y.dtype}")
    #    break

    return train_dataloader, test_dataloader


'''
Notes:
SGD with weight decay seems to slow down learning too much, don't use it
SGD with no LR result in jumping around towards the end
ADAM is very lr sensitive
Using pretrained weights results in getting stuck often
SGD with Nesterov momentum is worse
'''


def train_model(model, train_dataloader: DataLoader, test_dataloader: DataLoader, num_epochs: int, use_nesterov=False) -> bool:
    if num_epochs <= 0:
        logger.info("num epoch <= 0, skipping training")
        return False

    loss_fn = nn.CrossEntropyLoss(weight=train_dataloader.dataset.get_weights().to(device))

    '''
    # this paper shows SGD is hard to beat:
    # Closing the Generalization Gap of Adaptive Gradient Methods in Training Deep Neural Networks
    # They also did GRID search to optimise the parameters
    optimizer = torch.optim.SGD(model.parameters(), lr=0.1, momentum=0.9) # best so far
    optimizer = torch.optim.SGD(model.parameters(), lr=0.05, momentum=0.9) # best so far
    '''

    '''
    # use parameters from this paper: High-speed hyperparameter optimization for deep ResNet modelsin image recognition (2021)
    # but weight decay is bad
    #optimizer = torch.optim.SGD(model.parameters(), lr=0.1, momentum=0.9, weight_decay=5e-4)

    # parameters from this paper: Efficient ResNets: Residual Network Design (2023)
    # this seems to be the best paper that tested different optimisers, lr schedulers and along with the parameters
    # also see https://github.com/Nikunj-Gupta/Efficient_ResNets/blob/master/resnet_configs/config.yaml
    # with some changes such as no weight decay, sa it seems to worsen training for me
    optimizer = torch.optim.SGD(model.parameters(), lr=0.1, momentum=0.9)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=num_epochs, eta_min=optimizer.param_groups[0]["lr"] * 0.1)
    '''

    if use_nesterov:
        # nesterov with cosine annealing gives the best model, but trains slower
        optimizer = torch.optim.SGD(model.parameters(), lr=0.01, momentum=0.9, nesterov=True)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=num_epochs,
                                                               eta_min=optimizer.param_groups[0]["lr"] * 0.2)
    else:
        # adamw is best
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4, weight_decay=2.5e-2)

    train_stats_list = []
    test_stats_list = []

    logger.info(f"start training, num epochs={num_epochs}")

    # print the optimizer params
    optimizer_params = {k: v for (k, v) in optimizer.param_groups[0].items() if isinstance(v, (int, float))}
    logger.info(f"optimizer={optimizer.__class__.__name__}, params={optimizer_params}")

    # print the lr scheduler params
    try:
        logger.info(f"lr scheduler={scheduler.__class__.__name__}, states={scheduler.state_dict()}")
    except NameError:
        pass

    start = time.time()

    for e in range(num_epochs):
        epoch = e + 1
        train_epoch_stats = EpochStats()
        test_epoch_stats = EpochStats()
        train_stats_list.append(train_epoch_stats)
        test_stats_list.append(test_epoch_stats)

        # write log every 100 epoch
        should_log = epoch % 100 == 0

        # if should_log:
        logger.info(f"-----------------------------------------"
                    f" Epoch {epoch:>4d}  [lr={optimizer.param_groups[0]['lr']:.6f}] "
                    f"-----------------------------------------")

        train(train_dataloader, model, loss_fn, optimizer, train_epoch_stats, epoch, should_log)
        test(test_dataloader, model, loss_fn, test_epoch_stats, epoch, should_log)

        try:
            scheduler.step()
        except NameError:
            pass

    elapsed_min = (time.time() - start) / 60
    hours, minutes = divmod(elapsed_min, 60)
    logger.info(f"training complete! time taken: {hours:.0f}h {minutes:.0f}m")

    # create a df for the epoch data
    epoch_df = pd.DataFrame()
    epoch_df["epoch"] = [i + 1 for i in range(num_epochs)]
    epoch_df["trainLoss"] = [s.loss for s in train_stats_list]
    epoch_df["testLoss"] = [s.loss for s in test_stats_list]
    epoch_df["trainCount"] = [s.count for s in train_stats_list]
    epoch_df["trainNumCorrect"] = [s.num_correct for s in train_stats_list]
    epoch_df["testCount"] = [s.count for s in test_stats_list]
    epoch_df["testNumCorrect"] = [s.num_correct for s in test_stats_list]

    epoch_df.to_csv('epoch.tsv', sep='\t', index=False)

    epoch_df["testAccuracy"] = (epoch_df["testNumCorrect"] / epoch_df["testCount"]) * 100

    fig, ax_array = plt.subplots(2, 1, squeeze=False)
    fig.set_size_inches(10, 8)
    loss_ax = ax_array[0][0]
    acc_ax = ax_array[1][0]

    epoch_df.plot(ax=loss_ax, x="epoch", y="trainLoss", label="train loss", color="deepskyblue")
    epoch_df.plot(ax=loss_ax, x="epoch", y="testLoss", label="test loss", color="violet")

    epoch_df.plot(ax=acc_ax, x="epoch", y="testAccuracy", label="test accuracy%", color="red", style="--")

    loss_ax.set_xlim(0, num_epochs)
    loss_ax.set_ylim(0, 1)
    loss_ax.set_ylabel("loss")
    acc_ax.set_xlim(0, num_epochs)
    acc_ax.set_ylabel("accuracy %")
    plt.savefig('epoch.png', bbox_inches='tight')

    return True


def train(dataloader, model, loss_fn, optimizer, epoch_stats, epoch, should_log):
    size = len(dataloader.dataset)
    model.train()
    for batch, (image_x, linear_x, y) in enumerate(dataloader):
        image_x = image_x.to(device)
        linear_x = linear_x.to(device)
        y = y.to(device)
        # Compute prediction error
        pred = model(image_x, linear_x)

        # set weight to ignore low base qual bases
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        epoch_stats.update_loss(loss.item(), y.size(dim=0))

        # NOTE: need to add softmax
        pred = torch.nn.functional.softmax(pred, dim=1)
        epoch_stats.update_accuracies(pred, y)

        '''
        if should_log and ((batch + 1) % 10 == 0):
            loss = loss.item()
            current = (batch + 1) * y.size(dim=0)
            logger.debug(f"batch {batch} loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
        '''

    epoch_stats.log("Train", epoch)


def test(dataloader, model, loss_fn, epoch_stats, epoch, should_log):
    # size = len(dataloader.dataset)
    # num_batches = len(dataloader)
    model.eval()

    with torch.no_grad():
        for batch, (image_x, linear_x, y) in enumerate(dataloader):
            image_x = image_x.to(device)
            linear_x = linear_x.to(device)
            y = y.to(device)
            pred = model(image_x, linear_x)
            # pred = torch.full(y.shape, 1.0)

            # set weight to ignore low base qual bases
            # loss_fn.weight = float(1.0)

            epoch_stats.update_loss(loss_fn(pred, y).item(), y.size(dim=0))

            # NOTE: need to add softmax
            pred = torch.nn.functional.softmax(pred, dim=1)
            epoch_stats.update_accuracies(pred, y)
            # logger.debug(f"test size: x.size(0)")
            # logger.debug(f"true error: {t_err}, false error: {f_err}")

    epoch_stats.log("Test ", epoch)


def draw_confusion_matrix(dataloader, model, output_path='confusion_matrix.png'):
    from sklearn.metrics import ConfusionMatrixDisplay

    # After epoch stats logging, compute predictions for confusion matrix
    true_indices, pred_indices = [], []
    model.eval()
    with torch.no_grad():
        for image_x, linear_x, y in dataloader:
            image_x = image_x.to(device)
            linear_x = linear_x.to(device)
            y = y.to(device)
            pred = model(image_x, linear_x)
            pred_indices.extend(torch.argmax(torch.nn.functional.softmax(pred, dim=1), dim=1).cpu().numpy())
            true_indices.extend(torch.argmax(y, dim=1).cpu().numpy())

    fig, ax = plt.subplots(figsize=(10, 10))
    ConfusionMatrixDisplay.from_predictions(
        [cancer_types[i] for i in true_indices],
        [cancer_types[i] for i in pred_indices],
        # display_labels=cancer_types,
        cmap=plt.cm.Blues,
        xticks_rotation="vertical",
        ax=ax
    )
    plt.title('Confusion Matrix')
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Confusion matrix saved at {output_path}")


def append_dropout(model, rate):
    for name, module in model.named_children():
        if len(list(module.children())) > 0:
            append_dropout(module, rate)
        if isinstance(module, nn.ReLU):
            new = nn.Sequential(module, nn.Dropout2d(p=rate))
            setattr(model, name, new)


#
def init_model(cnn_dropout_rate, linear_dropout_rate):
    # load the model
    # 5
    return vcuppa_model.VCuppaModel(cnn_dropout_rate, linear_dropout_rate)


def train_main(sample_tsv, purple_root, driver_tsv, snv_tsv, epochs, batch_size, cnn_dropout_rate,
               linear_dropout_rate, test_fraction, use_nesterov, starting_model):
    df = pd.read_csv(sample_tsv, sep="\t")
    driver_df = pd.read_csv(driver_tsv, sep="\t")
    snv_df = pd.read_csv(snv_tsv, sep="\t")
    #snv_df = snv_df[snv_df["Category"] != "gen_pos"]
    snv_df = snv_df[snv_df["Category"] == "sig"]

    df = filter_df(df)

    #df = df.head(50).reset_index(drop=True)

    # remove any sample that is not driver df
    df = df[df["sampleId"].isin(driver_df.columns)].reset_index(drop=True)
    df = df[df["sampleId"].isin(snv_df.columns)].reset_index(drop=True)

    df["circosPngPath"] = purple_root + "/" + df["sampleId"] + "/" + df["sampleId"] + ".circos.png"

    # add all the driver to the df
    df["drivers"] = df["sampleId"].apply(lambda x: driver_df[x].to_numpy())
    df["snvs"] = df["sampleId"].apply(lambda x: snv_df[x].to_numpy())

    if starting_model:
        logger.info(f"starting model: {starting_model}")
        model = torch.jit.load(starting_model, map_location=torch.device('cpu'))
    else:
        model = init_model(cnn_dropout_rate, linear_dropout_rate)

    model = model.to(device)

    train_dataloader, test_dataloader = create_dataloader(df, batch_size=batch_size,
                                                          augment=True,
                                                          test_fraction=test_fraction)

    train_model(model, train_dataloader, test_dataloader, num_epochs=int(epochs), use_nesterov=use_nesterov)

    # Generate confusion matrix
    draw_confusion_matrix(test_dataloader, model)

    # save the model, always save in cpu
    model = model.to("cpu")
    torch.save(model.state_dict(), 'model.pth')

    # also use torch script to save it once more
    # load using
    # model = torch.jit.load('model_scripted.pt')
    model_scripted = torch.jit.script(model)  # Export to TorchScript
    model_scripted.save('model_scripted.pt')  # Save


def main():
    import argparse
    parser = argparse.ArgumentParser(description="train hrd predictor")
    parser.add_argument('--sample_tsv', help='input tsv file', required=True)
    parser.add_argument('--purple_root', help='path to purple plots', required=True)
    parser.add_argument('--driver_tsv', help='path to drivers tsv', required=True)
    parser.add_argument('--snv_tsv', help='path to snv tsv', required=True)
    parser.add_argument('--epochs', help='number epochs', type=int, default=200)
    parser.add_argument('--batch_size', help='batch size', type=int, default=25)
    parser.add_argument('--no_augment', help='disable augment training set by random rotation', action='store_true')
    parser.add_argument('--cnn_dropout_rate', help='dropout rate for the CNN', type=float, default=0.1)
    parser.add_argument('--linear_dropout_rate', help='dropout rate for the linear part', type=float, default=0.4)
    parser.add_argument('--test_fraction', help='amount of data used for testing', type=float, default=0.2)
    parser.add_argument('--use_nesterov', help='use SGD with nesterov instead of adamW', action='store_true')
    parser.add_argument('--starting_model', help='starting from this model instead of make a new one', default=None)
    args = parser.parse_args()

    logger.info(f"using {device} device, image size={IMAGE_SIZE}")
    logger.info(f"training cnn hrd, epochs={args.epochs}, batch_size={args.batch_size}, " +
                f"cnn_dropout_rate={args.cnn_dropout_rate}, linear_dropout_rate={args.linear_dropout_rate}, "
                f"test_fraction={args.test_fraction}, use_nesterov={args.use_nesterov}")
    logger.info(f"starting_model={args.starting_model}")
    logger.info(f"sample_tsv={args.sample_tsv}, driver_tsv={args.driver_tsv}, snv_tsv={args.snv_tsv}")
    logger.info(f"purple_root={args.purple_root}")
    train_main(args.sample_tsv, args.purple_root, args.driver_tsv, args.snv_tsv, args.epochs, args.batch_size,
               args.cnn_dropout_rate, args.linear_dropout_rate,
               args.test_fraction, args.use_nesterov, args.starting_model)


if __name__ == "__main__":
    main()
