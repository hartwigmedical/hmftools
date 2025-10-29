package com.hartwig.hmftools.geneutils.mapping;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

// Creates the index "image" file for the GATK BWA-MEM JNI wrapper, from a FASTA file.
// This image file is required for aligning with the GATK BWA-MEM JNI wrapper.
// It is basically a packed representation of the files you get from running `bwa index`
public class CreateGatkBwaMemIndex
{
    private static final String CFG_INPUT_FASTA = "input_fasta";
    private static final String CFG_OUTPUT_FILE = "output_file";

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        configBuilder.addPath(CFG_INPUT_FASTA, true, "Path to input FASTA file");
        configBuilder.addConfigItem(CFG_OUTPUT_FILE, false, "Path to output file");
        configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);

        configBuilder.checkAndParseCommandLine(args);

        String inputFasta = configBuilder.getValue(CFG_INPUT_FASTA);
        String outputFile = configBuilder.getValue(CFG_OUTPUT_FILE);
        if(outputFile == null)
        {
            outputFile = inputFasta + ".img";
        }

        loadAlignerLibrary(configBuilder.getValue(BWA_LIB_PATH));

        BwaMemIndex.createIndexImageFromFastaFile(inputFasta, outputFile);
    }
}
