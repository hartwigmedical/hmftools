package com.hartwig.hmftools.vchord;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import static com.hartwig.hmftools.vchord.VChordConfig.registerConfig;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.vchord.ImmutableVChordPrediction;
import com.hartwig.hmftools.common.vchord.VChordPrediction;
import com.hartwig.hmftools.common.vchord.VChordPredictionFile;

import org.apache.logging.log4j.Level;
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

public class VChordApplication
{
    public static final Logger LOGGER = LogManager.getLogger(VChordApplication.class);

    private final VChordConfig mConfig;

    public VChordApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new VChordConfig(configBuilder);
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

        // purity
        final double purity = PurityContextFile.read(mConfig.getPurpleDir(), mConfig.getSampleId()).bestFit().purity();

        LOGGER.info("loaded circos image({}), purity({})", circosPngPath, purity);

        ImmutableVChordPrediction.Builder vChordPredBuilder = ImmutableVChordPrediction.builder();

        // try to load the pytorch model
        Path modelPath = Paths.get(mConfig.getModelPath());
        try(Model model = Model.newInstance(modelPath.getFileName().toString(), Device.cpu()))
        {
            model.load(modelPath.getParent());

            // we make one prediction per cancer type
            for(HrdCancerType cancerType : HrdCancerType.values())
            {
                VChordInput input = new VChordInput(circosPng, cancerType, purity);

                // perform prediction on the image
                try(Predictor<VChordInput, Classifications> predictor = model.newPredictor(new PurplePlotTranslater()))
                {
                    double pred = predictor.predict(input).item(0).getProbability();

                    LOGGER.printf(Level.INFO, "%s hrd score: %.3f", cancerType.name().toLowerCase(), pred);

                    switch(cancerType)
                    {
                        case BREAST:
                            vChordPredBuilder.breastCancerHrdScore(pred);
                            break;
                        case OVARIAN:
                            vChordPredBuilder.ovarianCancerHrdScore(pred);
                            break;
                        case PANCREATIC:
                            vChordPredBuilder.pancreaticCancerScore(pred);
                            break;
                        case PROSTATE:
                            vChordPredBuilder.prostateCancerScore(pred);
                            break;
                        case OTHER:
                            vChordPredBuilder.otherCancerScore(pred);
                            break;
                    }
                }
            }
        }
        catch(TranslateException e)
        {
            throw new RuntimeException(e);
        }

        VChordPrediction vChordPrediction = vChordPredBuilder.build();

        // write to output
        VChordPredictionFile.write(VChordPredictionFile.generateFilename(mConfig.getOutputDir(), mConfig.getSampleId()),
                mConfig.getSampleId(), vChordPrediction);

        LOGGER.info("vChord complete");
    }

    public static void main(final String[] args) throws MalformedModelException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("vChord");
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        VersionInfo versionInfo = new VersionInfo("vchord.version");
        LOGGER.info("build timestamp: {}", versionInfo.buildTime().format(ISO_ZONED_DATE_TIME));

        new VChordApplication(configBuilder).run();
    }
}
