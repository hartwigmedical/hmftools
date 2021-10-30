package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.svtools.germline.GermlineUtils.GM_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.variant.vcf.VCFHeader;

public class VcfUtils
{
    // VCF field identifiers - can these be sourced from hmf-common SV classes?
    public static final String QUAL = "QUAL";
    public static final String SR = "SR";
    public static final String BQ = "BQ";
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

    public static GenotypeIds parseVcfSampleIds(final VCFHeader header, final String referenceId, final String tumorId)
    {
        List<String> vcfSampleNames = header.getGenotypeSamples();

        int tumorOrdinal = -1;
        int referenceOrdinal = -1;
        String vcfTumorId = "";
        String vcfRefefenceId = "";

        for(int i = 0; i < vcfSampleNames.size(); ++i)
        {
            String vcfSampleName = vcfSampleNames.get(i);

            if(vcfSampleName.contains(tumorId))
            {
                vcfTumorId = vcfSampleNames.get(i);
                tumorOrdinal = i;
            }
            else if(!referenceId.isEmpty() && vcfSampleName.contains(referenceId))
            {
                vcfRefefenceId = vcfSampleNames.get(i);
                referenceOrdinal = i;
            }
        }

        if(tumorOrdinal < 0 || (!referenceId.isEmpty() && referenceOrdinal < 0))
        {
            GM_LOGGER.error("missing sample names in VCF: {}", vcfSampleNames);
            return null;
        }


        return new GenotypeIds(referenceOrdinal, tumorOrdinal, vcfRefefenceId, vcfTumorId);
    }


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
