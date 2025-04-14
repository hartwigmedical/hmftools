package com.hartwig.hmftools.geneutils.mapping;

import static com.hartwig.hmftools.common.bam.CigarUtils.calcCigarAlignedLength;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.mapping.SequenceTester.FLD_SEQ_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.umccr.java.hellbender.utils.bwa.BwaMemAligner;
import org.umccr.java.hellbender.utils.bwa.BwaMemAlignment;
import org.umccr.java.hellbender.utils.bwa.BwaMemIndex;

public class BwaSeqTester
{
    private final SeqTestConfig mConfig;
    private final BwaMemAligner mAligner;

    public BwaSeqTester(final SeqTestConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;

        String refGenomeImageFile = configBuilder.getValue(REF_GENOME) + ".img";

        mAligner = initialiseBwaAligner(refGenomeImageFile);
    }

    private BwaMemAligner initialiseBwaAligner(final String refGenomeImageFile)
    {
        if(refGenomeImageFile.isEmpty() || !Files.exists(Paths.get(refGenomeImageFile)))
            System.exit(1);

        try
        {
            BwaMemIndex index = new BwaMemIndex(refGenomeImageFile);
            BwaMemAligner aligner = new BwaMemAligner(index);

            // log default values()

            GU_LOGGER.debug("BWA align options: ");
            GU_LOGGER.debug("  BWA Bandwidth: {}", aligner.getBandwidthOption());
            GU_LOGGER.debug("  BWA MatchScore: {}", aligner.getMatchScoreOption());
            GU_LOGGER.debug("  BWA OutputScoreThreshold: {}", aligner.getOutputScoreThresholdOption());
            GU_LOGGER.debug("  BWA SplitFactor: {}", aligner.getSplitFactorOption());
            GU_LOGGER.debug("  BWA Mismatch: {}", aligner.getXADropRatio());
            GU_LOGGER.debug("  BWA IGapOpen: {}", aligner.getIGapOpenPenaltyOption());
            GU_LOGGER.debug("  BWA IGapOpenExtend: {}", aligner.getIGapExtendPenaltyOption());
            GU_LOGGER.debug("  BWA DGapOpen: {}", aligner.getDGapOpenPenaltyOption());
            GU_LOGGER.debug("  BWA DGapOpenExtend: {}", aligner.getDGapExtendPenaltyOption());
            GU_LOGGER.debug("  BWA DropRatio: {}", aligner.getDropRatioOption());
            GU_LOGGER.debug("  BWA ZDrop: {}", aligner.getZDropOption());
            GU_LOGGER.debug("  BWA XADropRatio: {}", aligner.getXADropRatio());
            GU_LOGGER.debug("  BWA MaxXAHits: {}", aligner.getMaxXAHitsOption());
            GU_LOGGER.debug("  BWA MaxXAHitsAlt: {}", aligner.getMaxXAHitsAltOption());

            aligner.setBandwidthOption(31); // matches Esvee / Sage indel calling
            // aligner.setMismatchPenaltyOption(5);

            /*
            aligner.setOutputScoreThresholdOption(0);
            aligner.setZDropOption(1);
            aligner.setXADropRatio((float)0.1);
            aligner.setMatchScoreOption(1); // default 1
            aligner.setMaxXAHitsOption(100); // default 5
            aligner.setMaxXAHitsAltOption(200); // default 200
            aligner.setSplitFactorOption((float)0.1); // default 1.5
            */

            return aligner;
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            return null;
        }
    }

    public void run(final List<SequenceInfo> sequenceInfos)
    {
        long startTimeMs = System.currentTimeMillis();

        try
        {
            String outputFile = mConfig.OutputDir + mConfig.OutputId + ".bwa_alignments.tsv";

            BufferedWriter writer = initialiseBwaWriter(outputFile);

            GU_LOGGER.info("running BWA alignment for {} sequences", sequenceInfos.size());

            for(SequenceInfo sequenceInfo : sequenceInfos)
            {
                List<List<BwaMemAlignment>> bwaAlignments = mAligner.alignSeqs(List.of(sequenceInfo.Sequence.getBytes()));

                List<BwaMemAlignment> bwaSeqAlignments = bwaAlignments.get(0);

                GU_LOGGER.trace("seqId({}) region({}) found {} alignments",
                        sequenceInfo.Id, sequenceInfo.Region, bwaSeqAlignments.size());

                writeAlignments(sequenceInfo, bwaSeqAlignments, writer);
            }

            writer.close();
        }
        catch(Exception e)
        {
            GU_LOGGER.error("error reading / writing: {}", e.toString());
            System.exit(1);
        }

        GU_LOGGER.info("BWA aligner complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private BufferedWriter initialiseBwaWriter(final String filename) throws IOException
    {
        BufferedWriter writer = createBufferedWriter(filename);

        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(FLD_SEQ_ID);

        sj.add("RefLocation");
        sj.add("MapQual");
        sj.add("Cigar");
        sj.add("AlignedBases");
        sj.add("Score");
        sj.add("Flags");
        sj.add("NMatches");
        sj.add("XaTag");
        sj.add("MdTag");

        SequenceInfo.addSequenceHeader(sj);

        writer.write(sj.toString());
        writer.newLine();

        return writer;
    }

    private void writeAlignments(
            final SequenceInfo sequenceInfo, final List<BwaMemAlignment> alignments, final BufferedWriter writer) throws IOException
    {
        for(BwaMemAlignment alignment : alignments)
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(String.valueOf(sequenceInfo.Id));

            int chrIndex = alignment.getRefId();
            String chromosome;

            if(chrIndex >= 0 && chrIndex < HumanChromosome.values().length)
            {
                chromosome = mConfig.RefGenomeVersion.versionedChromosome(HumanChromosome.values()[chrIndex].toString());
            }
            else
            {
                chromosome = String.valueOf(chrIndex);
            }

            ChrBaseRegion refLocation = new ChrBaseRegion(chromosome, alignment.getRefStart() + 1, alignment.getRefEnd());

            String cigar = alignment.getCigar();
            int alignedBases = calcCigarAlignedLength(cigar);

            sj.add(refLocation.toString());
            sj.add(String.valueOf(alignment.getMapQual()));
            sj.add(cigar);
            sj.add(String.valueOf(alignedBases));
            sj.add(String.valueOf(alignment.getAlignerScore()));
            sj.add(String.valueOf(alignment.getSamFlag()));
            sj.add(String.valueOf(alignment.getNMismatches()));
            sj.add(alignment.getXATag());
            sj.add(alignment.getMDTag());

            sequenceInfo.addSequenceInfo(sj);

            writer.write(sj.toString());
            writer.newLine();
        }
    }
}
