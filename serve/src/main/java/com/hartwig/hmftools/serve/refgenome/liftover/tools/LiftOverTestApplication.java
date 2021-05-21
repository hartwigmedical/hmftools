package com.hartwig.hmftools.serve.refgenome.liftover.tools;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.serve.ServeConfig;
import com.hartwig.hmftools.serve.ServeLocalConfigProvider;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

public class LiftOverTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(LiftOverTestApplication.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        ServeConfig config = ServeLocalConfigProvider.create();

        // X:1314965 is a hotspot in 37 on CRLF2 that can't be lifted over to 38.
        Interval interval1 = new Interval("chrX", 1314965, 1314965);
        from37To38(config, interval1);

        // This ATM variant has changed ref
        Interval interval2 = new Interval("chr11", 108312440, 108312440);
        from38To37(config, interval2);
    }

    private static void from37To38(@NotNull ServeConfig config, @NotNull Interval original) {
        LOGGER.debug("Starting interval 37 is {}", original);

        LiftOver liftOver37To38 = new LiftOver(new File(config.refGenome37To38Chain()));
        Interval result = liftOver37To38.liftOver(original);
        LOGGER.debug("Interval lifted from 37 to 38 is {}", result);

        if (result != null) {
            LiftOver liftOver38to37 = new LiftOver(new File(config.refGenome38To37Chain()));
            Interval backToOriginal = liftOver38to37.liftOver(result);
            LOGGER.debug("Interval lifted back to 37 is {}", backToOriginal);
        }
    }

    private static void from38To37(@NotNull ServeConfig config, @NotNull Interval original) {
        LOGGER.debug("Starting interval 38 is {}", original);

        LiftOver liftOver38to37 = new LiftOver(new File(config.refGenome38To37Chain()));
        Interval result = liftOver38to37.liftOver(original);
        LOGGER.debug("Interval lifted from 38 to 37 is {}", result);

        if (result != null) {
            LiftOver liftOver37To38 = new LiftOver(new File(config.refGenome37To38Chain()));
            Interval backToOriginal = liftOver37To38.liftOver(result);
            LOGGER.debug("Interval lifted back to 38 is {}", backToOriginal);
        }
    }
}
