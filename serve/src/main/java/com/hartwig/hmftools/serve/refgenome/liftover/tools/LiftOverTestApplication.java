package com.hartwig.hmftools.serve.refgenome.liftover.tools;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.serve.ServeConfig;
import com.hartwig.hmftools.serve.ServeLocalConfigProvider;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

public class LiftOverTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(LiftOverTestApplication.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        ServeConfig config = ServeLocalConfigProvider.create();

        // X:1314965 is a hotspot in 37 on CRLF2 that can't be lifted over to 38.
        Interval original = new Interval("chrX", 1314965, 1314965);
        LOGGER.debug("Starting interval 37 is {}", original);

        LiftOver liftOver37To38 = new LiftOver(new File(config.refGenome37To38Chain()));
        Interval result = liftOver37To38.liftOver(original);
        LOGGER.debug("Interval lifted from 37 to 38 is {}", result);

        LiftOver liftOver38to37 = new LiftOver(new File(config.refGenome38To37Chain()));
        if (result != null) {
            Interval backToOriginal = liftOver38to37.liftOver(result);
            LOGGER.debug("Interval lifted back to 37 is {}", backToOriginal);
        }
    }
}
