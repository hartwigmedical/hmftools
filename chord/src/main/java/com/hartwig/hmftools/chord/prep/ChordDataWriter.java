package com.hartwig.hmftools.chord.prep;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

import org.jetbrains.annotations.Nullable;

public class ChordDataWriter
{
    private final String mOutputFile;
    private final BufferedWriter mWriter;

    private static final String FLD_SAMPLE_ID = "sample_id";

    public ChordDataWriter(String outputFile) throws IOException
    {
        mOutputFile = outputFile;
        mWriter = initializeWriter();
    }

    private BufferedWriter initializeWriter() throws IOException
    {
        CHORD_LOGGER.info("Writing mutation context counts to: {}", mOutputFile);
        return FileWriterUtils.createBufferedWriter(mOutputFile, false);
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
