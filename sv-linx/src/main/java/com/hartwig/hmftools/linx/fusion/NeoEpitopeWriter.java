package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.switchStream;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.linx.types.LinkedPair;

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
    }

    public void initialiseSample(final String sampleId)
    {
        mSampleId = sampleId;

        // clear any cache
        mFusions.clear();
    }

    public void processFusionCandidate(
            final List<GeneAnnotation> breakendGenes1, final List<GeneAnnotation> breakendGenes2,
            final List<LinkedPair> traversedPairs, final DisruptionFinder disruptionFinder)
    {
        if(breakendGenes1.isEmpty() || breakendGenes2.isEmpty())
            return;

        // avoid writing duplicates of the same junction
        if(isDuplicate(breakendGenes1.get(0), breakendGenes2.get(0)))
            return;

        boolean fusionAdded = false;

        for (final GeneAnnotation gene1 : breakendGenes1)
        {
            for (final GeneAnnotation gene2 : breakendGenes2)
            {
                // can allow upstream to upstream in reverse order??
                if (gene1.isUpstream() == gene2.isUpstream())
                    continue;

                final GeneAnnotation upGene = gene1.isUpstream() ? gene1 : gene2;

                if(!isCandidateUpstreamGene(upGene))
                    continue;

                final GeneAnnotation downGene = upGene == gene1 ? gene2 : gene1;

                if(!hasValidTraversal(upGene, downGene, traversedPairs, gene1.isUpstream(), disruptionFinder))
                    continue;

                double avgJcn = (upGene.jcn() + downGene.jcn()) * 0.5;

                NeoEpitopeFusion fusion = new NeoEpitopeFusion(
                        upGene.StableId, upGene.GeneName, upGene.chromosome(), upGene.position(), upGene.orientation(), upGene.id(),
                        downGene.StableId, downGene.GeneName, downGene.chromosome(), downGene.position(), downGene.orientation(), downGene.id(),
                        avgJcn, upGene.insertSequence());

                writeData(fusion);

                if(!fusionAdded)
                {
                    mFusions.add(fusion);
                    fusionAdded = true;
                }
            }
        }
    }

    private boolean hasValidTraversal(
            final GeneAnnotation upGene, final GeneAnnotation downGene, final List<LinkedPair> traversedPairs,
            boolean fusionLowerToUpper, final DisruptionFinder disruptionFinder)
    {
        if(traversedPairs == null || traversedPairs.isEmpty())
            return true;

        int upGeneStrand = upGene.Strand;
        boolean isPrecodingUpstream = false;

        for(LinkedPair pair : traversedPairs)
        {
            // if going lower to upper, if the orientation of the first breakend in the pair is opposite to the strand of
            // the upstream gene, then the fusion direction for that pair is the same as a the upstream gene
            // otherwise it needs to be switched
            int fusionDirection = 0;

            if(fusionLowerToUpper)
            {
                fusionDirection = pair.firstBreakend().orientation() != upGeneStrand ? upGeneStrand : -upGeneStrand;
            }
            else
            {
                fusionDirection = pair.secondBreakend().orientation() != upGeneStrand ? upGeneStrand : -upGeneStrand;
            }

            // any invalid traversal causes this fusion to be entirely skipped from further analysis
            if(disruptionFinder.pairTraversesGene(pair, fusionDirection, isPrecodingUpstream))
                return false;
        }

        return true;
    }

    private boolean isCandidateUpstreamGene(final GeneAnnotation gene)
    {
        // must have at least 1 coding transcript
        return gene.transcripts().stream().anyMatch(x -> x.isCoding() && x.isDisruptive());
    }

    private boolean isDuplicate(final GeneAnnotation gene1, final GeneAnnotation gene2)
    {
        for(final NeoEpitopeFusion fusion : mFusions)
        {
            for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
            {
                if(fusion.SvIds[fs] == gene1.id() && fusion.SvIds[switchStream(fs)] == gene2.id()
                && fusion.Chromosomes[fs].equals(gene1.chromosome()) && fusion.Chromosomes[switchStream(fs)].equals(gene2.chromosome())
                && fusion.Positions[fs] == gene1.position() && fusion.Positions[switchStream(fs)] == gene2.position())
                {
                    return true;
                }
            }
        }

        return false;
    }

    private void writeData(final NeoEpitopeFusion fusion)
    {
        try
        {
            if(mFileWriter == null)
            {
                if(mIsMultiSample)
                {
                    String outputFileName = mOutputDir + "LNX_NEO_EPITOPES.csv";

                    mFileWriter = createBufferedWriter(outputFileName, false);

                    mFileWriter.write("sampleId");
                    mFileWriter.write(NeoEpitopeFusion.header());
                    mFileWriter.newLine();
                }
                else
                {
                    String outputFileName = NeoEpitopeFusion.generateFilename(mOutputDir, mSampleId);

                    mFileWriter = createBufferedWriter(outputFileName, false);
                    mFileWriter.write(NeoEpitopeFusion.header());
                    mFileWriter.newLine();
                }
            }

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
