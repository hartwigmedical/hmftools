package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_DEPTH_END;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_DEPTH_START;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.util.StringJoiner;

public final class CanonicalSpliceJunctionFile
{
    public static final String CANONICAL_SJ_FILE_ID = "canonical_splice_junc.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + CANONICAL_SJ_FILE_ID;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(FLD_GENE_ID)
                .add(FLD_GENE_NAME)
                .add(FLD_CHROMOSOME)
                .add(FLD_ALT_SJ_POS_START)
                .add(FLD_ALT_SJ_POS_END)
                .add(FLD_ALT_SJ_TYPE)
                .add(FLD_FRAG_COUNT)
                .add(FLD_DEPTH_START)
                .add(FLD_DEPTH_END)
                .add("TranscriptNames")
                .toString();
    }
}
