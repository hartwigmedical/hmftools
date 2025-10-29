package com.hartwig.hmftools.vcuppa;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.VersionInfo;
import com.hartwig.hmftools.common.vcuppa.ImmutableVCuppaPrediction;
import com.hartwig.hmftools.common.vcuppa.VCuppaPrediction;
import com.hartwig.hmftools.common.vcuppa.VCuppaPredictionFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import ai.djl.Device;
import ai.djl.MalformedModelException;
import ai.djl.Model;
import ai.djl.inference.Predictor;
import ai.djl.modality.Classifications;
import ai.djl.modality.cv.Image;
import ai.djl.modality.cv.ImageFactory;
import ai.djl.translate.TranslateException;

public class VCuppaApplication
{
    public static final Logger LOGGER = LogManager.getLogger(VCuppaApplication.class);

    private final VCuppaConfig mConfig;

    public VCuppaApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new VCuppaConfig(configBuilder);
    }

    public void run() throws MalformedModelException, IOException
    {
        Path circosPngPath = Path.of(String.format("%s/plot/%s.circos.png", mConfig.getPurpleDir(), mConfig.getSampleId()));
        Image circosPng;

        try
        {
            circosPng = ImageFactory.getInstance().fromFile(circosPngPath);
        }
        catch(IOException e)
        {
            LOGGER.error("unable to load circos png({}), exception({})", circosPngPath, e.toString());
            throw e;
        }

        // linear features
        LinearFeatures linearFeatures = DriversAndSignaturesFile.load(mConfig.getCuppaFeatureDefinitionsPath(),
                mConfig.getCuppaFeaturesPath());

        LOGGER.info("loaded circos image({}), linear features({})", circosPngPath, mConfig.getCuppaFeaturesPath());

        VCuppaInput input = new VCuppaInput(circosPng, linearFeatures);

        ImmutableVCuppaPrediction.Builder vCuppaPredBuilder = ImmutableVCuppaPrediction.builder();

        // try to load the pytorch model
        LOGGER.info("loading model from path({})", mConfig.getModelPath());
        Path modelPath = Paths.get(mConfig.getModelPath()).toAbsolutePath();
        try(Model model = Model.newInstance(modelPath.getFileName().toString(), Device.cpu()))
        {
            model.load(modelPath.getParent());

            // perform prediction on the image
            try(Predictor<VCuppaInput, Classifications> predictor = model.newPredictor(new VCuppaInputTranslater()))
            {
                Classifications classifications = predictor.predict(input);

                // add all items to the pred
                classifications.items().forEach(o -> vCuppaPredBuilder.addCancerTypePredictions(
                        ImmutableVCuppaPrediction.CancerTypePrediction.builder()
                                .cancerType(o.getClassName())
                                .probability(o.getProbability())
                                .build()
                ));

                Classifications.Classification best = classifications.best();
                LOGGER.info("{} predicted cancer type({}), prob({})", mConfig.getSampleId(), best.getClassName(), best.getProbability());
            }
        }
        catch(TranslateException e)
        {
            throw new RuntimeException(e);
        }

        VCuppaPrediction vCuppaPrediction = vCuppaPredBuilder.build();

        // write to output
        VCuppaPredictionFile.write(VCuppaPredictionFile.generateFilename(mConfig.getOutputDir(), mConfig.getSampleId()),
                mConfig.getSampleId(), vCuppaPrediction);

        LOGGER.info("vChord complete");
    }

    public static void main(final String[] args) throws MalformedModelException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("vCuppa");
        VCuppaConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        VersionInfo versionInfo = new VersionInfo("vcuppa.version");
        LOGGER.info("build timestamp: {}", versionInfo.buildTime().format(ISO_ZONED_DATE_TIME));

        new VCuppaApplication(configBuilder).run();
    }
}
