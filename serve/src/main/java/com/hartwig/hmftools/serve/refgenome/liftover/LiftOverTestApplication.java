package com.hartwig.hmftools.serve.refgenome.liftover;

import java.io.File;
import java.net.InetAddress;
import java.net.UnknownHostException;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

public class LiftOverTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(LiftOverTestApplication.class);

    public static void main(String[] args) throws UnknownHostException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String liftOverDir;
        if (hostname.toLowerCase().contains("datastore")) {
            liftOverDir = "/data/common/refgenomes/liftover";
        } else {
            liftOverDir = System.getProperty("user.home") + "/hmf/refgenomes/liftover";
        }

        String chain37To38 = liftOverDir + "/hg19ToHg38.over.chain";
        String chain38To37 = liftOverDir + "/hg38ToHg19.over.chain";

        Interval original = new Interval("chr17", 41197710, 41197710);
        LOGGER.debug("Starting interval is {}", original);

        LiftOver liftOver37To38 = new LiftOver(new File(chain37To38));
        Interval result = liftOver37To38.liftOver(original);
        LOGGER.debug("Interval lifted from 37 to 38 is {}", result);

        LiftOver liftOver38to37 = new LiftOver(new File(chain38To37));
        Interval backToOriginal = liftOver38to37.liftOver(result);
        LOGGER.debug("Interval lifted back to 38 is {}", backToOriginal);
    }
}
