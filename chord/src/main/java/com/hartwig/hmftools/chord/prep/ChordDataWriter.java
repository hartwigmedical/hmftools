package com.hartwig.hmftools.chord.prep;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

import org.jetbrains.annotations.Nullable;

public class ChordDataWriter
{
    public final ChordConfig mConfig;

    private final BufferedWriter mWriter;

    public static final String MUTATION_CONTEXTS_FILE_SUFFIX = ".chord.mutation_contexts.tsv";
    public static final String COHORT_FILE_PREFIX = "cohort";
    private static final String FLD_SAMPLE_ID = "sample_id";

    public ChordDataWriter(ChordConfig config) throws IOException
    {
        mConfig = config;

        mWriter = initializeWriter();
    }

    public static String formOutputFile(String outputDir, String sampleId, @Nullable String outputId)
    {
        String outputFile = outputDir + "/" + sampleId;

        if(outputId != null)
            outputFile += "." + outputId;

        outputFile += MUTATION_CONTEXTS_FILE_SUFFIX;

        return outputFile;
    }

    @VisibleForTesting
    public static String formOutputFile(ChordConfig config, String sampleId)
    {
        return formOutputFile(config.OutputDir, sampleId, config.OutputId);
    }

    private BufferedWriter initializeWriter() throws IOException
    {
        String sampleIdOrCohort = mConfig.isSingleSample() ?
                mConfig.SampleIds.get(0) :
                COHORT_FILE_PREFIX;

        String outputFile = formOutputFile(mConfig.OutputDir, sampleIdOrCohort, mConfig.OutputId);

        CHORD_LOGGER.info("Writing mutation context counts to: {}", outputFile);

        return FileWriterUtils.createBufferedWriter(outputFile, false);
    }

    public void writeHeader(List<MutContextCount> contextCountsFirstSample) throws IOException
    {
        List<String> headerStrings = new ArrayList<>();

        headerStrings.add(FLD_SAMPLE_ID);

        for(MutContextCount mutContext : contextCountsFirstSample)
            headerStrings.add(mutContext.mName);

        String line = String.join(TSV_DELIM, headerStrings);
        mWriter.write(line);
        mWriter.newLine();
    }

    public void writeValues(String sampleId, List<MutContextCount> contextCounts) throws IOException
    {
        List<String> lineStrings = new ArrayList<>();

        lineStrings.add(sampleId);

        for(MutContextCount mutContext : contextCounts)
            lineStrings.add(String.valueOf(mutContext.mCount));

        String line = String.join(TSV_DELIM, lineStrings);
        mWriter.write(line);
        mWriter.newLine();
    }

    public void close() throws IOException
    {
        mWriter.close();
    }
}
