package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.svtools.germline.GermlineUtils.GM_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;

public class VcfUtils
{
    // VCF field identifiers
    public static final String QUAL = "QUAL";
    public static final String SR = "SR";
    public static final String SRQ = "SRQ";
    public static final String VF = "VF";
    public static final String RP = "RP";
    public static final String RPQ = "RPQ";
    public static final String REF = "REF";
    public static final String BEID = "BEID";
    public static final String BEIDL = "BEIDL";
    public static final String HOMSEQ = "HOMSEQ";

    public static final String AS = "AS";
    public static final String CAS = "CAS";
    public static final String RAS = "RAS";

    public static final String BASRP = "BASRP";
    public static final String ASRP = "ASRP";
    public static final String SB = "SB";
    public static final String BVF = "BVF";
    public static final String REFPAIR = "REFPAIR";


    public static List<String> loadVcfFiles(final String vcfsFile)
    {
        List<String> vcfFiles = Lists.newArrayList();

        if (!Files.exists(Paths.get(vcfsFile)))
            return vcfFiles;

        try
        {
            vcfFiles = Files.readAllLines(new File(vcfsFile).toPath());

            GM_LOGGER.info("loaded {} VCF filenames", vcfFiles.size());
        }
        catch(IOException e)
        {
            GM_LOGGER.error("failed to load gene panel file({}): {}", vcfsFile, e.toString());
        }

        return vcfFiles;
    }
}
