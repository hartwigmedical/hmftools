import pandas as pd

import numpy as np
import math
import os, sys, time

import torch
from torch import nn
import torch.utils.data as data_utils

import torchvision
# use v2 as it claims to be faster
from torchvision.transforms import v2
from torchvision.io import read_image

import matplotlib.pyplot as plt

import logging

import hrd_model

logger = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s %(levelname)5s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logger.setLevel(logging.DEBUG)

IMAGE_SIZE = 512
NUM_CANCER_TYPES = 5

# select clinical.sampleId, clinical.primaryTumorLocation, chord.BRCA1, chord.BRCA2, chord.hrd, chord.hrStatus, chord.hrdType, chord.remarksHrStatus, chord.remarksHrdType from clinical, chord where clinical.sampleId = chord.sampleId and (clinical.primaryTumorLocation = 'Breast' or clinical.primaryTumorLocation = 'Ovary' or clinical.primaryTumorLocation = 'Fallopian tube' or clinical.primaryTumorLocation = "Prostate" or clinical.primaryTumorLocation = "Pancreas");

# Get cpu or gpu device for training.
device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"

# convert a df to tensor to be used in pytorch
def df_to_tensor(df):
    return torch.from_numpy(df.values).float().to(device)


def filter_df(df):
    # filter out stuff we do not need, we filter out anything that seem abiguous to Chord
    df = df[(df["hrStatus"] == "HR_DEFICIENT") | (df["hrStatus"] == "HR_PROFICIENT")]
    return df


# https://stackoverflow.com/questions/44865023/how-can-i-create-a-circular-mask-for-a-numpy-array
def create_circular_mask(h, w, center=None, radius=None):
    if center is None:  # use the middle of the image
        center = (int(w / 2), int(h / 2))
    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w - center[0], h - center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)

    mask = dist_from_center >= radius
    return torch.from_numpy(mask.astype(np.uint8))  # note we must cast to uint8


# create a transform object to convert image to pytorch tensor
def make_transform(image_size=IMAGE_SIZE):
    x = [v2.Resize(int(image_size * 1.1)),
          v2.CenterCrop(image_size),
          v2.ToDtype(torch.float32, scale=True)]
    return v2.Compose(x)

def image_to_tensor(input_png_path, circos_png_path, image_size=IMAGE_SIZE, transform=None):
    if transform is None:
        transform = make_transform(image_size)

    # input png we use gray scale image
    #input_png = transform(v2.functional.to_grayscale(read_image(input_png_path)))
    #input_png = transform(read_image(input_png_path))
    circos_png = transform(read_image(circos_png_path))
    #return torch.cat((input_png, circos_png), dim=0)
    return circos_png

def cancer_type_to_tensor(cancer_type, purity):
    a = [0.0, 0.0, 0.0, 0.0, 0.0]
    i = -1
    if isinstance(cancer_type, str):
        cancer_type = cancer_type.lower()
        if cancer_type == "breast":
            i = 0
        elif cancer_type == "ovary" or cancer_type == "ovarian":
            i = 1
        elif cancer_type == "pancreas" or cancer_type == "pancreatic":
            i = 2
        elif cancer_type == "prostate":
            i = 3
    if i != -1:
        a[i] = 1.0
    a[4] = purity
    #logger.info(f"cancer type: {cancer_type}, encoded: {a}")
    return torch.tensor(a, dtype=torch.float32)

class HrdDataset(data_utils.Dataset):

    def __init__(self, df, image_size, augment, hrd_sample_dup):
        start = time.time()

        # duplicate all certain HRD samples to enrich the size
        df = pd.concat([df] + [df[df["hrd"] >= 0.8]] * hrd_sample_dup, ignore_index=True)

        is_hrd = df["hrd"].astype(np.float32)

        num_hrd = len(df[df["hrStatus"] == "HR_DEFICIENT"])

        self.targets = [torch.tensor([x]) for x in is_hrd]
        self.image_tensors = []
        self.type_tensors = []

        transform = make_transform(image_size)

        logger.info(f"start loading images, augment={augment}")
        # load all the images into an array
        for idx, row in df.iterrows():
            self.image_tensors.append(image_to_tensor(row["inputPngPath"], row["circosPngPath"], image_size, transform))
            self.type_tensors.append(cancer_type_to_tensor(row["primaryTumorLocation"], row["purity"]))

        # rotation is done when the tensor is retrieved, this provides better
        # augmentation
        if augment:
            self.transform = v2.RandomRotation(degrees=180, fill=1.0)
        else:
            self.transform = None

        elapsed_sec = int(time.time() - start)
        minute, second = divmod(elapsed_sec, 60)
        logger.info(f"loading dataset of size {len(self)}, hrd={num_hrd}, took {minute:.0f}m {second:.0f}s")

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
        type_tensor = self.type_tensors[idx]
        target_tensor = self.targets[idx]
        return image_tensor, type_tensor, target_tensor


class EpochStats:
    def __init__(self):
        self.loss_sum = 0.0
        self.loss_count = 0
        self.true_pos = 0
        self.num_pos = 0
        self.true_neg = 0
        self.num_neg = 0

    @property
    def loss(self):
        return self.loss_sum / self.loss_count

    @property
    def accuracy(self):
        return (self.true_pos + self.true_neg) / (self.num_pos + self.num_neg)

    @property
    def correct(self):
        return self.true_pos + self.true_neg

    @property
    def count(self):
        return self.num_pos + self.num_neg

    @property
    def true_pos_rate(self):
        return self.true_pos / self.num_pos

    @property
    def true_neg_rate(self):
        return self.true_neg / self.num_neg

    def update_loss(self, loss_val, count):
        self.loss_sum += loss_val * count
        self.loss_count += count

    def update_accuracies(self, pred: torch.Tensor, target: torch.Tensor):
        # logger.debug(f"calc_error: pred: {pred}, target: {target}")

        # force prediction to 0 and 1, this is required after sigmoid
        pred = (pred > 0.5).int()
        target = (target > 0.5).int()
        self.num_pos += target.sum().item()
        self.num_neg += (1 - target).sum().item()

        self.true_pos += (pred * target).sum().item()

        # 1 - target gives the mask required to zero out all true target ones
        self.true_neg += ((1 - pred) * (1 - target)).sum().item()

    def log(self, name, epoch):
        logger.info(f"[{name} {epoch:>4d}]    loss: {self.loss:>8.5f}    "
                    f"acc:{self.accuracy * 100:>6.2f}% [{self.correct:>4d}/{self.count:>4d}]    "
                    f"TP:{self.true_pos_rate * 100:>6.2f}% [{self.true_pos:>4d}/{self.num_pos:>4d}]    "
                    f"TN:{self.true_neg_rate * 100:>6.2f}% [{self.true_neg:>4d}/{self.num_neg:>4d}]")


# from the input pandas dataframe create the train / test dataloaders requied for
# pytorch training.
# return (train_dataloader, test_dataloader)
def create_dataloader(df, image_size, batch_size, augment, hrd_sample_dup):
    # Using Skicit-learn to split data into training and testing sets
    from sklearn.model_selection import train_test_split

    # Split the data into training and testing sets
    train_df, test_df = train_test_split(df, test_size=0.25, random_state=None)

    # write out the train and test set
    train_df[["sampleId"]].to_csv("train_set.tsv.gz", sep="\t", index=False)
    test_df[["sampleId"]].to_csv("test_set.tsv.gz", sep="\t", index=False)

    logger.info(f"train size: {len(train_df)}, test size: {len(test_df)}, batch size: {batch_size}, augment: {augment}")

    # Create data loaders.
    train_dataset = HrdDataset(train_df, image_size, augment=augment, hrd_sample_dup=hrd_sample_dup)
    train_dataloader = data_utils.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

    # we do not augment testing data
    test_dataset = HrdDataset(test_df, image_size, augment=False, hrd_sample_dup=hrd_sample_dup)
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
def train_model(model, train_dataloader, test_dataloader, num_epochs, use_nesterov=False):
    # loss_fn = nn.BCELoss()

    # use this cause the model does not have a logit
    loss_fn = nn.BCEWithLogitsLoss()

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
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=num_epochs, eta_min=optimizer.param_groups[0]["lr"] * 0.2)
    else:
        # adamw is best
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4, weight_decay=2.5e-2)

    train_stats_list = []
    test_stats_list = []

    logger.info(f"start training, num epochs={num_epochs}")

    # print the optimizer params
    optimizer_params = {k:v for (k,v) in optimizer.param_groups[0].items() if isinstance(v, (int, float))}
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
    epoch_df["trainNumPos"] = [s.num_pos for s in train_stats_list]
    epoch_df["trainTruePos"] = [s.true_pos for s in train_stats_list]
    epoch_df["trainNumNeg"] = [s.num_neg for s in train_stats_list]
    epoch_df["trainTrueNeg"] = [s.true_neg for s in train_stats_list]
    epoch_df["testNumPos"] = [s.num_pos for s in test_stats_list]
    epoch_df["testTruePos"] = [s.true_pos for s in test_stats_list]
    epoch_df["testNumNeg"] = [s.num_neg for s in test_stats_list]
    epoch_df["testTrueNeg"] = [s.true_neg for s in test_stats_list]

    epoch_df.to_csv('epoch.tsv', sep='\t', index=False)

    epoch_df["testTP"] = epoch_df["testTruePos"] / epoch_df["testNumPos"] * 100
    epoch_df["testTN"] = epoch_df["testTrueNeg"] / epoch_df["testNumNeg"] * 100
    epoch_df["testAccuracy"] = (epoch_df["testTruePos"] + epoch_df["testTrueNeg"]) / (
                epoch_df["testNumPos"] + epoch_df["testNumNeg"]) * 100

    fig, ax_array = plt.subplots(2, 1, squeeze=False)
    fig.set_size_inches(10, 8)
    loss_ax = ax_array[0][0]
    acc_ax = ax_array[1][0]

    epoch_df.plot(ax=loss_ax, x="epoch", y="trainLoss", label="train loss", color="deepskyblue")
    epoch_df.plot(ax=loss_ax, x="epoch", y="testLoss", label="test loss", color="violet")

    epoch_df.plot(ax=acc_ax, x="epoch", y="testAccuracy", label="test accuracy%", color="red", style="--")
    epoch_df.plot(ax=acc_ax, x="epoch", y="testTP", label="test TP%", color="teal")
    epoch_df.plot(ax=acc_ax, x="epoch", y="testTN", label="test TN%", color="orange")

    loss_ax.set_xlim(0, num_epochs)
    loss_ax.set_ylim(0, 1)
    loss_ax.set_ylabel("loss")
    acc_ax.set_xlim(0, num_epochs)
    acc_ax.set_ylabel("accuracy %")

    plt.savefig('epoch.png', bbox_inches='tight')

    # make sure we have a good model
    final_test_acc = epoch_df["testAccuracy"][-1]
    final_test_tp = epoch_df["testTP"][-1]
    final_test_tn = epoch_df["testTN"][-1]

    if final_test_tp < final_test_tn:
        logger.warn(f"final testTP[{final_test_tp:>6.2f}] < testTN[{final_test_tn:>6.2f}], please retrain")
        return False
    return True


def train(dataloader, model, loss_fn, optimizer, epoch_stats, epoch, should_log):
    size = len(dataloader.dataset)
    model.train()
    for batch, (x1, x2, y) in enumerate(dataloader):
        x1 = x1.to(device)
        x2 = x2.to(device)
        y = y.to(device)
        # Compute prediction error
        pred = model(x1, x2)

        # set weight to ignore low base qual bases
        # loss_fn.weight = float(1.0)
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        epoch_stats.update_loss(loss.item(), y.size(dim=0))

        # sigmoid is needed since we use loss with logit
        sigmoid_pred = torch.sigmoid(pred)
        epoch_stats.update_accuracies(sigmoid_pred, y)

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
        for x1, x2, y in dataloader:
            x1 = x1.to(device)
            x2 = x2.to(device)
            y = y.to(device)
            pred = model(x1, x2)
            # pred = torch.full(y.shape, 1.0)

            # set weight to ignore low base qual bases
            # loss_fn.weight = float(1.0)

            epoch_stats.update_loss(loss_fn(pred, y).item(), y.size(dim=0))

            # NOTE: need to add sigmoid since we use loss with logit
            sigmoid_pred = torch.sigmoid(pred)
            epoch_stats.update_accuracies(sigmoid_pred, y)
            # logger.debug(f"test size: x.size(0)")
            # logger.debug(f"true error: {t_err}, false error: {f_err}")

    epoch_stats.log("Test ", epoch)


def append_dropout(model, rate):
    for name, module in model.named_children():
        if len(list(module.children())) > 0:
            append_dropout(module, rate)
        if isinstance(module, nn.ReLU):
            new = nn.Sequential(module, nn.Dropout2d(p=rate))
            setattr(model, name, new)

#
def init_model(dropout_rate):
    # load the model
    # 5
    return hrd_model.HrdModel(dropout_rate, NUM_CANCER_TYPES)

def train_main(sample_tsv, purple_root, epochs, batch_size, augment, dropout_rate, hrd_sample_dup, use_nesterov, starting_model):
    df = pd.read_csv(sample_tsv, sep="\t")
    df["inputPngPath"] = purple_root + "/" + df["sampleId"] + "/" + df["sampleId"] + ".input.png"
    df["circosPngPath"] = purple_root + "/" + df["sampleId"] + "/" + df["sampleId"] + ".circos.png"

    # load the purity
    df["purity"] = [pd.read_csv(f'{purple_root}/{s}/{s}.purple.purity.tsv', sep='\t')["purity"].iloc[0] for s in df["sampleId"]]

    df = filter_df(df)
    #df = df.head(5).reset_index(drop=True)

    if starting_model:
        logger.info(f"starting model: {starting_model}")
        model = torch.jit.load(starting_model, map_location=torch.device('cpu'))
    else:
        model = init_model(dropout_rate)

    model = model.to(device)

    train_dataloader, test_dataloader = create_dataloader(df, image_size=IMAGE_SIZE, batch_size=batch_size,
                                                          augment=augment, hrd_sample_dup=hrd_sample_dup)
    train_model(model, train_dataloader, test_dataloader, num_epochs=epochs, use_nesterov=use_nesterov)

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
    parser.add_argument('--epochs', help='number epochs', type=int, default=200)
    parser.add_argument('--batch_size', help='batch size', type=int, default=25)
    parser.add_argument('--no_augment', help='disable augment training set by random rotation', action='store_true')
    parser.add_argument('--dropout_rate', help='dropout rate', type=float, default=0.2)
    parser.add_argument('--hrd_sample_duplication', help='number to duplicate HRD samples by', type=int, default=7)
    parser.add_argument('--use_nesterov', help='use SGD with nesterov instead of adamW', action='store_true')
    parser.add_argument('--starting_model', help='starting from this model instead of make a new one', default=None)
    args = parser.parse_args()

    logger.info(f"using {device} device")
    logger.info(f"training cnn hrd, epochs={args.epochs}, batch_size={args.batch_size}, augment={not args.no_augment}, " +
          f"dropout_rate={args.dropout_rate}, hrd_sample_dup={args.hrd_sample_duplication}, use_nesterov={args.use_nesterov}, " +
          f"starting_model={args.starting_model}")
    train_main(args.sample_tsv, args.purple_root, args.epochs, args.batch_size, not args.no_augment, args.dropout_rate,
               args.hrd_sample_duplication, args.use_nesterov, args.starting_model)


if __name__ == "__main__":
    main()
