package com.hartwig.hmftools.geneutils.mapping;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.calcCigarAlignedLength;
import static com.hartwig.hmftools.common.blastn.BlastnMatch.calcSumBitScore;
import static com.hartwig.hmftools.common.blastn.BlastnRunner.BLAST_DB;
import static com.hartwig.hmftools.common.blastn.BlastnRunner.BLAST_TOOL;
import static com.hartwig.hmftools.common.blastn.BlastnRunner.registerBlastn;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.BLASTN_WORD_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.MIN_BLAST_ALIGNMENT_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.blastn.BlastnMatch;
import com.hartwig.hmftools.common.blastn.BlastnRunner;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.jetbrains.annotations.NotNull;

public class SequenceTester
{
    private static final String INPUT_FILE = "input_file";

    private final String mInputFile;
    private final String mOutputDir;
    private final String mOutputId;

    private final int mThreads;

    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenomeVersion;

    private final BwaMemAligner mAligner;
    private final BlastnRunner mBlastnRunner;

    private static final String FLD_SEQ_ID = "SeqId";
    private static final String FLD_SEQUENCE = "Sequence";

    public SequenceTester(final ConfigBuilder configBuilder)
    {
        String refGenomeFile = configBuilder.getValue(REF_GENOME);
        mRefGenome = loadRefGenome(refGenomeFile);
        mRefGenomeVersion = deriveRefGenomeVersion(mRefGenome);

        mInputFile = configBuilder.getValue(INPUT_FILE);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);
        mThreads = parseThreads(configBuilder);

        loadAlignerLibrary(configBuilder.getValue(BWA_LIB_PATH));

        String refGenomeImageFile = refGenomeFile + ".img";
        mAligner = initialiseBwaAligner(refGenomeImageFile);

        mBlastnRunner = initialiseBlastN(configBuilder.getValue(BLAST_TOOL), configBuilder.getValue(BLAST_DB));
    }
    
    private BwaMemAligner initialiseBwaAligner(final String refGenomeImageFile)
    {
        if(refGenomeImageFile.isEmpty() || !Files.exists(Paths.get(refGenomeImageFile)))
            System.exit(1);

        try
        {
            BwaMemIndex index = new BwaMemIndex(refGenomeImageFile);
            BwaMemAligner aligner = new BwaMemAligner(index);
            aligner.setBandwidthOption(31); // matches Esvee / Sage indel calling
            return aligner;
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            return null;
        }
    }

    private BlastnRunner initialiseBlastN(final String blastNTool, final String blastNDb)
    {
        if(blastNTool == null || blastNDb == null)
            return null;

        return new BlastnRunner.Builder()
                .withTask("megablast")
                .withPrefix(mOutputId)
                .withBlastDir(blastNTool)
                .withBlastDb(blastNDb)
                .withOutputDir(mOutputDir)
                .withKeepOutput(true)
                .withWordSize(BLASTN_WORD_SIZE)
                .withSubjectBestHit(true)
                .withNumThreads(mThreads)
                .build();
    }

    private class SequenceInfo
    {
        public final int Id;
        public final String Sequence;
        public final ChrBaseRegion Region;

        public SequenceInfo(final int id, final String sequence, final ChrBaseRegion region)
        {
            Id = id;
            Sequence = sequence;
            Region = region;
        }

        public String toString() { return format("id(%d) region(%ss) seq(%d:%s)", Id, Region, Sequence.length(), Sequence); }
    }

    public void run()
    {
        if(mInputFile == null || !Files.exists(Paths.get(mInputFile)))
        {
            GU_LOGGER.error("invalid input file");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        List<SequenceInfo> sequenceInfos = loadSequenceInfo();

        runBwaAligner(sequenceInfos);

        runBlastN(sequenceInfos);

        GU_LOGGER.info("sequence tester complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<SequenceInfo> loadSequenceInfo()
    {
        List<SequenceInfo> sequenceInfos = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(mInputFile));
            lines.remove(0); // Skip header

            int chrIndex = 0;
            int posStartIndex = 1;
            int posEndIndex = 2;
            int sequenceIndex = 3;

            int seqId = 0;
            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM);

                String chromosome = values[chrIndex];
                String startString = values[posStartIndex];
                String endString = values[posEndIndex];
                String sequence = values.length > sequenceIndex ? values[sequenceIndex] : null;

                ChrBaseRegion region = null;
                if(chromosome.length() > 0 && startString.length() > 0 && endString.length() > 0)
                {
                    region = new ChrBaseRegion(chromosome, Integer.parseInt(startString), Integer.parseInt(endString));
                }

                // Get sequence from region if sequence is empty
                if(sequence == null || sequence.length() == 0)
                {
                    sequence = mRefGenome.getBaseString(chromosome, region.start(), region.end());
                }

                sequenceInfos.add(new SequenceInfo(seqId, sequence, region));
                seqId++;
            }

            GU_LOGGER.info("loaded {} sequences from file({})", sequenceInfos.size(), mInputFile);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error reading input file({}): {}", mInputFile, e.toString());
            System.exit(1);
        }

        return sequenceInfos;
    }

    private void runBwaAligner(final List<SequenceInfo> sequenceInfos)
    {
        long startTimeMs = System.currentTimeMillis();

        try
        {
            String outputFile = mOutputDir + mOutputId + ".bwa_alignments.tsv";

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
        sj.add(FLD_CHROMOSOME).add(FLD_POS_START).add(FLD_POS_END);
        sj.add(FLD_SEQUENCE);

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
                chromosome = mRefGenomeVersion.versionedChromosome(HumanChromosome.values()[chrIndex].toString());
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

            if(sequenceInfo != null)
            {
                sj.add(sequenceInfo.Region.Chromosome);
                sj.add(String.valueOf(sequenceInfo.Region.start()));
                sj.add(String.valueOf(sequenceInfo.Region.end()));
            }
            else
            {
                sj.add(null);
                sj.add(null);
                sj.add(null);
            }

            sj.add(sequenceInfo.Sequence);

            writer.write(sj.toString());
            writer.newLine();
        }
    }

    private void runBlastN(final List<SequenceInfo> sequenceInfos)
    {
        //long startTimeMs = System.currentTimeMillis();

        try
        {
            String outputFile = mOutputDir + mOutputId + ".blastn_results.tsv";

            BufferedWriter writer = initialiseBlastNWriter(outputFile);

            GU_LOGGER.info("running BlastN for {} sequences", sequenceInfos.size());

            Map<Integer,String> sequencesMap = Maps.newHashMap();

            for(SequenceInfo sequenceInfo : sequenceInfos)
            {
                sequencesMap.put(sequenceInfo.Id, sequenceInfo.Sequence);
            }

            Multimap<Integer,BlastnMatch> sequenceResults = mBlastnRunner.run(sequencesMap);

            for(SequenceInfo sequenceInfo : sequenceInfos)
            {
                Collection<BlastnMatch> matches = sequenceResults.get(sequenceInfo.Id);
                writeBlastNMatches(sequenceInfo, matches, writer);
            }

            writer.close();
        }
        catch(Exception e)
        {
            GU_LOGGER.error("error calling BlastN: {}", e.toString());
            System.exit(1);
        }

        // GU_LOGGER.info("BlastN complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private BufferedWriter initialiseBlastNWriter(final String filename) throws IOException
    {
        BufferedWriter writer = createBufferedWriter(filename);

        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(FLD_SEQ_ID);

        sj.add("ResultCount");
        sj.add("SumBitScore");
        sj.add("AlignedLength");
        sj.add("RefLocation");
        sj.add("BitScore");
        sj.add("Mismatchs");
        sj.add("GapOpen");

        sj.add(FLD_CHROMOSOME).add(FLD_POS_START).add(FLD_POS_END);
        sj.add(FLD_SEQUENCE);

        writer.write(sj.toString());
        writer.newLine();

        return writer;
    }

    private static int MAX_BLASTN_RESULTS = 20;

    private void writeBlastNMatches(
            final SequenceInfo sequenceInfo, final Collection<BlastnMatch> matches, final BufferedWriter writer) throws IOException
    {
        if(matches.isEmpty())
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(String.valueOf(sequenceInfo.Id));
            sj.add("0");
            sj.add("-1");
            sj.add("none");
            sj.add("-1");
            sj.add("-1");
            sj.add("-1");

            if(sequenceInfo != null)
            {
                sj.add(sequenceInfo.Region.Chromosome);
                sj.add(String.valueOf(sequenceInfo.Region.start()));
                sj.add(String.valueOf(sequenceInfo.Region.end()));
            }
            else
            {
                sj.add(null);
                sj.add(null);
                sj.add(null);
            }

            sj.add(sequenceInfo.Sequence);

            writer.write(sj.toString());
            writer.newLine();
            return;
        }

        double sumBitScore = calcSumBitScore(matches, MIN_BLAST_ALIGNMENT_LENGTH);
        int resultCount = matches.size();

        List<BlastnMatch> sortedMatches = matches.stream().collect(Collectors.toList());

        Collections.sort(sortedMatches, Comparator.comparingDouble(x -> -x.BitScore));

        for(int i = 0; i < min(sortedMatches.size(), MAX_BLASTN_RESULTS); ++i)
        {
            BlastnMatch match = sortedMatches.get(i);
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(String.valueOf(sequenceInfo.Id));

            sj.add(String.valueOf(resultCount));
            sj.add(String.valueOf(sumBitScore));
            sj.add(String.valueOf(match.AlignmentLength));

            String chromosome = "";
            String chrNameStr = "Homo sapiens chromosome ";
            int chrIndex = match.SubjectTitle.indexOf(chrNameStr);

            // chromosome 12, GRCh38
            // chromosome 6 genomic
            // chromosome 14 unlocalized

            if(chrIndex >= 0)
            {
                int chrEndIndex = match.SubjectTitle.indexOf(", GRCh");

                if(chrEndIndex == -1)
                    chrEndIndex = match.SubjectTitle.indexOf(" genomic");

                if(chrEndIndex == -1)
                    chrEndIndex = match.SubjectTitle.indexOf(" unlocalized");

                if(chrEndIndex >= 0)
                    chromosome = match.SubjectTitle.substring(chrNameStr.length() + chrIndex, chrEndIndex);
            }

            ChrBaseRegion refLocation = new ChrBaseRegion(chromosome, match.SubjectAlignStart, match.SubjectAlignEnd);

            sj.add(refLocation.toString());
            sj.add(String.valueOf(match.BitScore));
            sj.add(String.valueOf(match.NumMismatch));
            sj.add(String.valueOf(match.NumGapOpenings));

            if(sequenceInfo != null)
            {
                sj.add(sequenceInfo.Region.Chromosome);
                sj.add(String.valueOf(sequenceInfo.Region.start()));
                sj.add(String.valueOf(sequenceInfo.Region.end()));
            }
            else
            {
                sj.add(null);
                sj.add(null);
                sj.add(null);
            }

            sj.add(sequenceInfo.Sequence);

            writer.write(sj.toString());
            writer.newLine();
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(INPUT_FILE, true, "Input regions and sequences");

        configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);
        registerBlastn(configBuilder, false);

        addOutputOptions(configBuilder, false);
        ConfigUtils.addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addRefGenomeFile(configBuilder, true);

        configBuilder.checkAndParseCommandLine(args);

        SequenceTester sequenceTester = new SequenceTester(configBuilder);

        sequenceTester.run();
    }
}
