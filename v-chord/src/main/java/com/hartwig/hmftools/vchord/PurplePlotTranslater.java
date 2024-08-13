package com.hartwig.hmftools.vchord;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ai.djl.modality.Classifications;
import ai.djl.modality.cv.Image;
import ai.djl.modality.cv.transform.CenterCrop;
import ai.djl.modality.cv.transform.Resize;
import ai.djl.modality.cv.transform.ToTensor;
import ai.djl.modality.cv.util.NDImageUtils;
import ai.djl.ndarray.NDArray;
import ai.djl.ndarray.NDList;
import ai.djl.translate.Transform;
import ai.djl.translate.Translator;
import ai.djl.translate.TranslatorContext;

/*
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
    input_png = v2.functional.to_grayscale(transform(read_image(input_png_path)))
    circos_png = transform(read_image(circos_png_path))
    return torch.cat((input_png, circos_png), dim=0)
 */
public class PurplePlotTranslater implements Translator<VChordInput, Classifications>
{
    @Override
    public NDList processInput(TranslatorContext ctx, VChordInput input) {
        //
        // Convert Image to NDArray

        // input.toNDArray(ctx.getNDManager(), Image.Flag.GRAYSCALE);

        NDArray purpleCircosArray = imageTransform(ctx, input.purpleCircos);

        // see https://djl.ai/docs/pytorch/pytorch-djl-ndarray-cheatsheet.html
        NDArray cancerTypeArray = cancerTypeAndPurityToTensor(ctx, input.cancerType, input.purity);

        return new NDList(purpleCircosArray, cancerTypeArray);
    }

    @Override
    public Classifications processOutput(TranslatorContext ctx, NDList list)
    {
        // Create a Classifications with the output probabilities
        double hrdScore = sigmoid(list.singletonOrThrow().getFloat(0));
        return new Classifications(List.of("hrd"), List.of(hrdScore));
    }

    // perform the image transform
    NDArray imageTransform(TranslatorContext ctx, Image image)
    {
        List<Transform> transforms = new ArrayList<>();
        transforms.add(new Resize((int)(VChordConstants.IMAGE_SIZE * 1.1)));
        transforms.add(new CenterCrop(VChordConstants.IMAGE_SIZE, VChordConstants.IMAGE_SIZE));
        transforms.add(new ToTensor());

        NDArray imageArray = image.toNDArray(ctx.getNDManager()); // convert to NDArray

        for(Transform t : transforms)
        {
            imageArray = t.transform(imageArray);
        }
        return imageArray;
    }

    /*
    def cancer_type_to_tensor(cancer_type, purity):
        cancer_type = cancer_type.lower()
        if cancer_type == "breast":
            i = 0
        elif cancer_type == "ovary" or cancer_type == "ovarian":
            i = 1
        elif cancer_type == "pancreas" or cancer_type == "pancreatic":
            i = 2
        elif cancer_type == "prostate":
            i = 3
        else:
            i = 4
        a = [0.0, 0.0, 0.0, 0.0, 0.0]
        a[i] = 1.0
        a[4] = purity
        #logger.info(f"cancer type: {cancer_type}, encoded: {a}")
        return torch.tensor(a, dtype=torch.float32)
     */
    NDArray cancerTypeAndPurityToTensor(TranslatorContext ctx, HrdCancerType cancerType, double purity)
    {
        int i;
        switch(cancerType)
        {
            case BREAST:
                i = 0;
                break;
            case OVARIAN:
                i = 1;
                break;
            case PANCREATIC:
                i = 2;
                break;
            case PROSTATE:
                i = 3;
                break;
            default:
                i = 4;
                break;
        }
        float[] array = new float[5];
        array[i] = 1.0f;
        array[4] = (float)purity;
        return ctx.getNDManager().create(array);
    }

    // sigmoid function
    public static double sigmoid(float x)
    {
        // Numerically-stable sigmoid function
        if (x >= 0)
        {
            double z = Math.exp(-x);
            return 1 / (1 + z);
        } else
        {
            double z = Math.exp(x);
            return z / (1 + z);
        }
    }
}
