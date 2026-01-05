package com.hartwig.hmftools.sage.candidate;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TIER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.common.VariantReadContext;

public class CandidateWriter
{
    private final SageCallConfig mConfig;
    private final BufferedWriter mWriter;

    private static final String CANDIDATES_FILE_ID = ".candidates.tsv.gz";

    public CandidateWriter(final SageCallConfig config)
    {
        mConfig = config;
        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter); }

    private BufferedWriter initialiseWriter()
    {
        if(!SageCallConfig.LogCandidates)
            return null;

        String outputVcf = mConfig.Common.OutputFile;
        String fileName = outputVcf.replace(".vcf.gz", CANDIDATES_FILE_ID);

        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_CHROMOSOME).add(FLD_POSITION).add(FLD_REF).add(FLD_ALT).add(FLD_TIER);

            sj.add("CoreBases").add("CoreCigar");
            sj.add("FullMatches").add("CoreMatches").add("LowQualInCore");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise candidates writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeCandidates(final List<Candidate> candidates)
    {
        if(mWriter == null)
            return;

        try
        {
            for(Candidate candidate : candidates)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(candidate.chromosome());
                sj.add(String.valueOf(candidate.position()));
                sj.add(candidate.variant().ref());
                sj.add(candidate.variant().alt());
                sj.add(String.valueOf(candidate.tier()));

                sj.add(candidate.readContext().coreStr());
                sj.add(candidate.readContext().readCigar());

                sj.add(String.valueOf(candidate.fullMatchSupport()));
                sj.add(String.valueOf(candidate.coreMatchSupport()));

                mWriter.write(sj.toString());
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write candidate: {}", e.toString());
        }
    }

    private static String readContextInfo(final VariantReadContext readContext)
    {
        return format("%s-%s", readContext.coreStr(), readContext.readCigar());
    }
}
