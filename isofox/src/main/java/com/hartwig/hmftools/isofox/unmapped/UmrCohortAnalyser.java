package com.hartwig.hmftools.isofox.unmapped;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.EXPRESSION_SCOPE_GENE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.START_STR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.startEndStr;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.UNMAPPED_READS;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.UMR_NO_MATE;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.baseQualsFromString;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class UmrCohortAnalyser
{
    private final CohortConfig mConfig;

    // other config
    private final RnaExpressionMatrix mGeneExpression;

    // map of chromosomes to a map of splice-boundary keys to a map of samples to a list of unmapped reads
    private final Map<String,Map<String,Map<String,List<UnmappedRead>>>> mUnmappedReads;

    private final LineElementMatcher mLineElementMatcher;
    private final boolean mCombineFrequencies;
    private final boolean mGroupBySequence;

    private final BufferedWriter mWriter;
    private final BufferedWriter mBlatWriter;
    private final BlatMatcher mBlatMatcher;
    private int mSequenceId;

    private static final String GENE_EXPRESSION_FILE = "gene_expression_file";
    private static final String COMBINE_FREQUENCIES = "combine_frequencies";
    private static final String GROUP_BY_SEQUENCE = "group_by_sequence";
    private static final String SV_VCF_FILE = "sv_vcf";
    private static final String WRITE_BLAT_FILE = "write_blat";
    private static final String BLAT_RESULTS_FILE = "blat_results_file";

    private static final String COHORT_FILE_ID = "unmapped_reads.csv";
    private static final String COHORT_RESULTS_FILE_ID = "unmapped_blat_results.csv";
    private static final int MIN_BLAT_SEQEUNCE_LENGTH = 20;
    private static final int REF_EXONIC_BASE_LENGTH = 20;

    public UmrCohortAnalyser(final CohortConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mUnmappedReads = Maps.newHashMap();

        mGeneExpression = configBuilder.hasValue(GENE_EXPRESSION_FILE) ?
                new RnaExpressionMatrix(configBuilder.getValue(GENE_EXPRESSION_FILE), EXPRESSION_SCOPE_GENE) : null;

        mLineElementMatcher = configBuilder.hasValue(LINX_DIR_CFG) && configBuilder.hasValue(SV_VCF_FILE) ?
                new LineElementMatcher(configBuilder.getValue(LINX_DIR_CFG), configBuilder.getValue(SV_VCF_FILE)) : null;

        mCombineFrequencies = configBuilder.hasFlag(COMBINE_FREQUENCIES);
        mGroupBySequence = configBuilder.hasFlag(GROUP_BY_SEQUENCE);
        mBlatMatcher = configBuilder.hasValue(BLAT_RESULTS_FILE) ? new BlatMatcher(configBuilder.getValue(BLAT_RESULTS_FILE)) : null;

        mWriter = initialiseWriter();
        mBlatWriter = configBuilder.hasFlag(WRITE_BLAT_FILE) ? initialiseBlatWriter() : null;
        mSequenceId = 0;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(GENE_EXPRESSION_FILE, false, "Gene expression file for cohort");
        configBuilder.addFlag(COMBINE_FREQUENCIES, "Determine cohort frequencies for slice candidates");
        configBuilder.addFlag(GROUP_BY_SEQUENCE, "Form consensus soft-clip sequences");
        configBuilder.addFlag(WRITE_BLAT_FILE, "Write fasta file for BLAT consensue sequence search");
        configBuilder.addPath(SV_VCF_FILE, false, "Structural variant VCF");
        configBuilder.addConfigItem(BLAT_RESULTS_FILE, false, "BLAT results file");
    }

    public void processSampleFiles()
    {
        if(mBlatMatcher != null)
        {
            ISF_LOGGER.info("matching consensus sequences to BLAT results");
            processCohortFile();
            return;

        }

        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, UNMAPPED_READS, filenames))
            return;

        int nextLog = 100000;

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path umrFile = filenames.get(i);

            loadFile(sampleId, umrFile);

            if(mCombineFrequencies)
            {
                int totalUmrCount = mUnmappedReads.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();

                if(totalUmrCount >= nextLog)
                {
                    ISF_LOGGER.debug("cached unmapped-read count({})", totalUmrCount);
                    nextLog += 100000;
                }
            }
            else
            {
                if(mLineElementMatcher != null)
                    mLineElementMatcher.findMatches(sampleId, mUnmappedReads);

                writeUnmappedReads();
                mUnmappedReads.clear();
                mSequenceId = 0; // reset for each sample
            }

            if(i > 0 && (i % 100) == 0)
            {
                ISF_LOGGER.debug("processed {} samples", i);
            }
        }

        if(mCombineFrequencies)
        {
            int totalUmrCount = mUnmappedReads.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();
            ISF_LOGGER.info("loaded {} unmapped-read records", totalUmrCount);

            ISF_LOGGER.info("writing cohort unmapped-read");

            writeUnmappedReads();
        }

        closeBufferedWriter(mWriter);
        closeBufferedWriter(mBlatWriter);
    }

    private void loadFile(final String sampleId, final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);

            lines.remove(0);

            int geneIdIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int geneNameIndex = fieldsIndexMap.get(FLD_GENE_NAME);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int readIndex = fieldsIndexMap.get("ReadId");
            int posStartIndex = fieldsIndexMap.get("ReadStart");
            int posEndIndex = fieldsIndexMap.get("ReadEnd");
            int spliceTypeIndex = fieldsIndexMap.get("SpliceType");
            int scLengthIndex = fieldsIndexMap.get("SoftClipLength");
            int scSideIndex = fieldsIndexMap.get("SoftClipSide");
            int abqIndex = fieldsIndexMap.get("AvgBaseQual");
            int transIndex = fieldsIndexMap.get("TransName");
            int exonRankIndex = fieldsIndexMap.get("ExonRank");
            int exonBoundaryIndex = fieldsIndexMap.get("ExonBoundary");
            int exonDistIndex = fieldsIndexMap.get("ExonDistance");
            int scBasesIndex = fieldsIndexMap.get("SoftClipBases");
            int scBaseQualsIndex = fieldsIndexMap.get("SoftClipBaseQuals");
            Integer mateIndex = fieldsIndexMap.get("MateCoords");
            Integer matchesSuppIndex = fieldsIndexMap.get("MatchesChimeric");

            for(String data : lines)
            {
                final String[] values = data.split(DELIMITER);

                UnmappedRead umRead = new UnmappedRead(
                        values[readIndex],
                        new ChrBaseRegion(values[chrIndex], Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex])),
                        Integer.parseInt(values[scLengthIndex]), Integer.parseInt(values[scSideIndex]),
                        Double.parseDouble(values[abqIndex]), values[geneIdIndex], values[geneNameIndex], values[transIndex],
                        Integer.parseInt(values[exonRankIndex]), Integer.parseInt(values[exonBoundaryIndex]),
                        Integer.parseInt(values[exonDistIndex]), values[spliceTypeIndex], values[scBasesIndex],
                        baseQualsFromString(values[scBaseQualsIndex], values[scBasesIndex].length()),
                        values[mateIndex], Boolean.parseBoolean(values[matchesSuppIndex]));

                addUnmappedRead(sampleId, umRead);
            }

            ISF_LOGGER.debug("sample({}) loaded {} unmapped-read records", sampleId, lines.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load unmapped-read file({}): {}", filename.toString(), e.toString());
        }
    }

    private void addUnmappedRead(final String sampleId, final UnmappedRead umRead)
    {
        Map<String,Map<String,List<UnmappedRead>>> chrUmrs = mUnmappedReads.get(umRead.ReadRegion.Chromosome);

        if(chrUmrs == null)
        {
            chrUmrs = Maps.newHashMap();
            mUnmappedReads.put(umRead.ReadRegion.Chromosome, chrUmrs);
        }

        String umrKey = umRead.positionKey();
        Map<String,List<UnmappedRead>> umrKeyList = chrUmrs.get(umrKey);

        if(umrKeyList == null)
        {
            umrKeyList = Maps.newHashMap();
            chrUmrs.put(umrKey, umrKeyList);
        }

        List<UnmappedRead> sampleReads = umrKeyList.get(sampleId);

        if(sampleReads == null)
        {
            sampleReads = Lists.newArrayList();
            umrKeyList.put(sampleId, sampleReads);
        }

        sampleReads.add(umRead);
    }

    private BufferedWriter initialiseWriter()
    {
        final String outputFile = mBlatMatcher != null ?
                mConfig.formCohortFilename(COHORT_RESULTS_FILE_ID) : mConfig.formCohortFilename(COHORT_FILE_ID);

        ISF_LOGGER.info("writing cohort file {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            if(mBlatMatcher != null)
                return writer;

            writer.write("SampleId,FragmentCount,UnpairedCount,Chromosome,GeneName,TransName,ExonRank,SpliceType");
            writer.write(",ExonBoundary,ExonDistance,SoftClipSide,AvgBaseQual");
            writer.write(",SoftClipBases,GeneTPM,HasChimericMatch");

            if(mCombineFrequencies)
                writer.write(",CohortFrequency");

            if(mLineElementMatcher != null)
                writer.write(",SvLinxMatches");

            if(mGroupBySequence)
                writer.write(",SequenceId,ExactMatchReads,AssignedTotal,HighQualPerc");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to initialise cohort unmapped reads file({}): {}", outputFile, e.toString());
            return null;
        }
    }

    private BufferedWriter initialiseBlatWriter()
    {
        final String outputFile = mConfig.formCohortFilename("blat_sequences.fa");

        try
        {
            return createBufferedWriter(outputFile, false);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to initialise BLAT seqeuence file({}): {}", outputFile, e.toString());
            return null;
        }
    }

    private void writeUnmappedReads()
    {
        if(mWriter == null)
            return;

        for(Map.Entry<String,Map<String,Map<String,List<UnmappedRead>>>> chrEntry : mUnmappedReads.entrySet())
        {
            for(Map<String,List<UnmappedRead>> umrSampleMap : chrEntry.getValue().values())
            {
                int sampleCount = umrSampleMap.size();

                for(Map.Entry<String,List<UnmappedRead>> sampleEntry : umrSampleMap.entrySet())
                {
                    String sampleId = sampleEntry.getKey();
                    List<UnmappedRead> umReads = sampleEntry.getValue();

                    if(mGroupBySequence)
                        writeLocationSequences(sampleId, umReads);
                    else
                        writeLocationReads(sampleId, umReads, sampleCount);
                }
            }
        }
    }

    private void writeLocationReads(final String sampleId, final List<UnmappedRead> umReads, int cohortSampleCount)
    {
        try
        {
            UnmappedRead firstRead = umReads.get(0);
            String chromosome = firstRead.ReadRegion.Chromosome;

            String matchedSVsInfo = mLineElementMatcher != null ? mLineElementMatcher.formUmrMatchString(firstRead.positionKey()) : "";

            // de-dup reads by fragment and across genes sharing the same exon boundary using readId
            Set<String> readIds = Sets.newHashSet();
            Set<String> genes = Sets.newHashSet();

            double geneTpm = 0;
            boolean hasSuppMatch = false;
            int unpairedCount = 0;

            for(UnmappedRead umRead : umReads)
            {
                if(!genes.contains(umRead.GeneName))
                {
                    genes.add(umRead.GeneName);

                    if(mGeneExpression != null)
                        geneTpm += mGeneExpression.getExpression(umRead.GeneId, sampleId);
                }

                readIds.add(umRead.ReadId);
                hasSuppMatch |= umRead.MatchesChimeric;

                if(umRead.MateCoords.equals(UMR_NO_MATE))
                    ++unpairedCount;
            }

            StringJoiner genesStr = new StringJoiner(ITEM_DELIM);
            genes.forEach(x -> genesStr.add(x));

            int fragmentCount = readIds.size();

            double avgBaseQual = umReads.stream().mapToDouble(x -> x.AvgBaseQual).sum() / umReads.size();

            mWriter.write(String.format("%s,%d,%d,%s,%s,%s,%d,%s",
                    sampleId, fragmentCount, unpairedCount, chromosome, genesStr.toString(), firstRead.TransName,
                    firstRead.ExonRank, firstRead.SpliceType));

            mWriter.write(String.format(",%d,%d,%s,%.1f,%s,%4.3e,%s",
                    firstRead.ExonBoundary, firstRead.ExonDistance, startEndStr(firstRead.ScSide),
                    avgBaseQual, firstRead.ScBases, geneTpm, hasSuppMatch));

            if(mCombineFrequencies)
            {
                mWriter.write(String.format(",%d", cohortSampleCount));
            }

            if(mLineElementMatcher != null)
            {
                mWriter.write(String.format(",%s", matchedSVsInfo));
            }

            mWriter.newLine();

        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write cohort unmapped reads file: {}", e.toString());
        }
    }

    private void writeLocationSequences(final String sampleId, final List<UnmappedRead> umReads)
    {
        try
        {
            SequenceTracker sequenceTracker = new SequenceTracker();
            umReads.forEach(x -> sequenceTracker.processRead(x));
            sequenceTracker.reconcileSequences();

            UnmappedRead firstRead = umReads.get(0);

            String chromosome = firstRead.ReadRegion.Chromosome;

            String matchedSVsInfo = mLineElementMatcher != null ? mLineElementMatcher.formUmrMatchString(firstRead.positionKey()) : "";

            Set<String> genes = Sets.newHashSet();

            double geneTpm = 0;

            for(UnmappedRead umRead : umReads)
            {
                if(!genes.contains(umRead.GeneName))
                {
                    genes.add(umRead.GeneName);

                    if(mGeneExpression != null)
                        geneTpm += mGeneExpression.getExpression(umRead.GeneId, sampleId);
                }
            }

            StringJoiner uniqueGenes = new StringJoiner(ITEM_DELIM);
            genes.forEach(x -> uniqueGenes.add(x));
            String genesStr = uniqueGenes.toString();

            List<ConsensusSequence> sequences = sequenceTracker.getSequences();

            for(ConsensusSequence sequence : sequences)
            {
                final List<UnmappedRead> assignedReads = sequence.allMatchedReads();

                // de-dup reads by fragment and across genes sharing the same exon boundary using readId
                Set<String> readIds = Sets.newHashSet();
                boolean hasSuppMatch = false;
                int unpairedCount = 0;

                for(UnmappedRead umRead : assignedReads)
                {
                    readIds.add(umRead.ReadId);
                    hasSuppMatch |= umRead.MatchesChimeric;

                    if(umRead.MateCoords.equals(UMR_NO_MATE))
                        ++unpairedCount;
                }

                int sequenceId = mSequenceId++;

                int fragmentCount = readIds.size();

                double avgBaseQual = assignedReads.stream().mapToDouble(x -> x.AvgBaseQual).sum() / assignedReads.size();

                mWriter.write(String.format("%s,%d,%d,%s,%s,%s,%d,%s",
                        sampleId, fragmentCount, unpairedCount, chromosome, genesStr, firstRead.TransName,
                        firstRead.ExonRank, firstRead.SpliceType));

                String consensusString = sequence.getSequenceString();

                mWriter.write(String.format(",%d,%d,%s,%.1f,%s,%4.3e,%s",
                        firstRead.ExonBoundary, firstRead.ExonDistance, startEndStr(firstRead.ScSide),
                        avgBaseQual, consensusString, geneTpm, hasSuppMatch));

                if(mLineElementMatcher != null)
                {
                    mWriter.write(String.format(",%s", matchedSVsInfo));
                }

                mWriter.write(String.format(",%d,%d,%.2f,%.2f",
                        sequenceId, sequence.exactMatchReads().size(), sequence.assignedReadTotal(), sequence.calcHighQualPercent()));

                mWriter.newLine();

                if(mBlatWriter != null && consensusString.length() >= MIN_BLAT_SEQEUNCE_LENGTH)
                {
                    mBlatWriter.write(String.format(">%s_%d", sampleId, sequenceId));
                    mBlatWriter.newLine();
                    mBlatWriter.write(consensusString);
                    mBlatWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write cohort unmapped reads file: {}", e.toString());
        }
    }

    private void processCohortFile()
    {
        final String inputFile = mConfig.formCohortFilename(COHORT_FILE_ID);

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(inputFile));

            String line = fileReader.readLine();
            String[] columnHeaders = line.split(DELIMITER);

            for(String string : columnHeaders)
            {
                mWriter.write(String.format("%s,", string));
            }

            mWriter.write("ExonicBases," + BlatResult.csvHeader());
            mWriter.newLine();

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIMITER);

            int sampleIdIndex = fieldsIndexMap.get("SampleId");
            int seqIdIndex = fieldsIndexMap.get("SequenceId");
            int chrIndex = fieldsIndexMap.get("Chromosome");
            int exonBoundaryIndex = fieldsIndexMap.get("ExonBoundary");
            int scSideIndex = fieldsIndexMap.get("SoftClipSide");

            List<BlatResult> sampleBlatResults = null;
            String currentSample = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DELIMITER, -1);

                String sampleId = values[sampleIdIndex];

                if(!currentSample.equals(sampleId))
                {
                    if(sampleBlatResults != null && !sampleBlatResults.isEmpty())
                    {
                        ISF_LOGGER.error("sample({}) has {} unmatched BLAT results", currentSample, sampleBlatResults.size());
                    }

                    currentSample = sampleId;
                    sampleBlatResults = Lists.newArrayList();
                    List<BlatResult> blatResults = mBlatMatcher.getSampleBlatResults(sampleId);

                    if(blatResults != null)
                        sampleBlatResults.addAll(blatResults);
                }

                if(sampleBlatResults.isEmpty())
                    continue;

                int sequenceId = Integer.parseInt(values[seqIdIndex]);

                BlatResult blatResult = mBlatMatcher.findBlatSequenceMatch(sampleBlatResults, sequenceId);

                if(blatResult == null) // ignore sequences without a blat result
                    continue;

                sampleBlatResults.remove(blatResult);

                for(String string : values)
                {
                    mWriter.write(String.format("%s,", string));
                }

                String exonicBases = getExonicBases(
                        values[chrIndex], Integer.parseInt(values[exonBoundaryIndex]), values[scSideIndex]);

                mWriter.write(String.format("%s,%s", exonicBases, blatResult.toCsv()));
                mWriter.newLine();
            }

            mWriter.close();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load cohort unmapped reads file({}): {}", inputFile.toString(), e.toString());
        }
    }

    private String getExonicBases(final String chromosome, int exonBoundary, final String scSide)
    {
        if(mConfig.RefGenome == null)
            return "";

        boolean isStart = scSide.equals(START_STR);

        int posStart = isStart ? exonBoundary : exonBoundary - REF_EXONIC_BASE_LENGTH + 1;
        int posEnd = isStart ? exonBoundary + REF_EXONIC_BASE_LENGTH + 1 : exonBoundary;
        String refBases = mConfig.RefGenome.getBaseString(chromosome, posStart, posEnd);

        return isStart ? Nucleotides.reverseComplementBases(refBases) : refBases;
    }

}
