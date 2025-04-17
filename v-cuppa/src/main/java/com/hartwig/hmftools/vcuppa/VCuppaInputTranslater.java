package com.hartwig.hmftools.vcuppa;

import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import ai.djl.modality.Classifications;
import ai.djl.modality.cv.Image;
import ai.djl.modality.cv.transform.Resize;
import ai.djl.modality.cv.transform.ToTensor;
import ai.djl.ndarray.NDArray;
import ai.djl.ndarray.NDList;
import ai.djl.translate.Transform;
import ai.djl.translate.Translator;
import ai.djl.translate.TranslatorContext;

/*
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
 */

public class VCuppaInputTranslater implements Translator<VCuppaInput, Classifications>
{
    public static final Logger LOGGER = LogManager.getLogger(VCuppaInputTranslater.class);

    @Override
    public NDList processInput(TranslatorContext ctx, VCuppaInput input) {

        // Convert Image to NDArray
        NDArray purpleCircosArray = imageTransform(ctx, input.purpleCircos());

        // LOGGER.info("image tensor: {}", purpleCircosArray.get(1).get(200).toFloatArray());

        // see https://djl.ai/docs/pytorch/pytorch-djl-ndarray-cheatsheet.html
        NDArray linearFeatureArray = ctx.getNDManager().create(input.linearFeatures().toFloatArray());

        return new NDList(purpleCircosArray, linearFeatureArray);
    }

    @Override
    public Classifications processOutput(TranslatorContext ctx, NDList list)
    {
        // LOGGER.info("raw pred: {}", list.singletonOrThrow());

        // Create a Classifications with the output probabilities
        // we need to apply softmax to the NDList
        NDArray propabilities = list.singletonOrThrow().softmax(0);
        return new Classifications(CancerTypes.cancerTypes, propabilities);
    }

    // perform the image transform
    NDArray imageTransform(TranslatorContext ctx, Image image)
    {
        List<Transform> transforms = new ArrayList<>();
        transforms.add(new Resize(VCuppaConstants.IMAGE_SIZE));
        transforms.add(new ToTensor());

        NDArray imageArray = image.toNDArray(ctx.getNDManager()); // convert to NDArray

        for(Transform t : transforms)
        {
            imageArray = t.transform(imageArray);
        }
        return imageArray;
    }
}
