package com.hartwig.hmftools.ckb;

import java.io.IOException;
import java.net.InetAddress;

import com.hartwig.hmftools.ckb.datamodel.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.reader.CkbJsonReader;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class CkbImporterTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbImporterTestApp.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running CKB Importer on '{}'", hostname);

        String ckbDir;
        Integer maxFilesToReadPerType = 10;

        if (hostname.toLowerCase().contains("datastore")) {
            ckbDir = "/data/common/dbs/ckb/210129_flex_dump";
        } else {
            ckbDir = System.getProperty("user.home") + "/hmf/projects/serve/ckb";
        }

        CkbJsonDatabase ckbDatabase = CkbJsonReader.read(ckbDir, maxFilesToReadPerType);

        LOGGER.info("Complete!");
    }
}
