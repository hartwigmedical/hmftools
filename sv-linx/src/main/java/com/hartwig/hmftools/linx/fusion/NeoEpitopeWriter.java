package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.isIrrelevantSameGene;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.validFusionTranscript;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;

import org.apache.commons.compress.utils.Lists;

public class NeoEpitopeWriter
{
    private final String mOutputDir;

    private boolean mIsMultiSample;
    private BufferedWriter mFileWriter;
    private String mSampleId;
    private final List<NeoEpitopeFusion> mFusions;

    public NeoEpitopeWriter(final String outputDir, boolean isMultiSample)
    {
        mOutputDir = outputDir;
        mIsMultiSample = isMultiSample;
        mFileWriter = null;
        mSampleId = null;
        mFusions = Lists.newArrayList();

        if(mIsMultiSample)
            intialiseMultiSampleWriter();
    }

    public void initialiseSample(final String sampleId)
    {
        mSampleId = sampleId;

        // clear any cache
        mFusions.clear();

        if(!mIsMultiSample)
        {

        }

    }

    public void processFusionCandidate(final List<GeneAnnotation> breakendGenes1, final List<GeneAnnotation> breakendGenes2)
    {
        for (final GeneAnnotation gene1 : breakendGenes1)
        {
            boolean startUpstream = gene1.isUpstream();

            final Transcript trans1 = gene1.transcripts().stream().filter(Transcript::isCanonical).findFirst().orElse(null);

            if (trans1 == null)
                continue;

            for (final GeneAnnotation gene2 : breakendGenes2)
            {
                boolean endUpstream = gene2.isUpstream();

                // can allow upstream to upstream in reverse order??
                if (startUpstream == endUpstream)
                    continue;

                NeoEpitopeFusion fusion = new NeoEpitopeFusion(
                        gene1.StableId, gene1.GeneName, gene1.chromosome(), gene1.position(), gene1.orientation(), gene1.id(),
                        gene2.StableId, gene2.GeneName, gene2.chromosome(), gene2.position(), gene2.orientation(), gene2.id());

                // TO-DO - avoid writing duplicates of the same junction
                if(isDuplicate(fusion))
                    continue;

                mFusions.add(fusion);

                writeData(fusion);
            }
        }
    }

    private boolean isDuplicate(final NeoEpitopeFusion fusion)
    {
        // mFusions
        return false;
    }

    private void intialiseMultiSampleWriter()
    {
        if(mOutputDir.isEmpty())
            return;

        try
        {
            String outputFileName = mOutputDir + "LNX_NEO_EPITOPES.csv";

            mFileWriter = createBufferedWriter(outputFileName, false);

            mFileWriter.write("sampleId");
            mFileWriter.write(NeoEpitopeFusion.header());
            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing neo-epitope output file: {}", e.toString());
        }
    }

    private void writeData(final NeoEpitopeFusion fusion)
    {
        if(mFileWriter == null)
            return;

        try
        {
            if(mIsMultiSample)
            {
                mFileWriter.write(String.format("%s",mSampleId));
            }

            mFileWriter.write(NeoEpitopeFusion.toString(fusion));
            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing neo-epitope output file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

}
