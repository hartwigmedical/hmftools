package com.hartwig.hmftools.geneutils.mapping;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.blastn.BlastnMatch.calcSumBitScore;
import static com.hartwig.hmftools.common.blastn.BlastnMatch.isPrimaryBlastnMatch;
import static com.hartwig.hmftools.common.blastn.BlastnRunner.BLAST_DB;
import static com.hartwig.hmftools.common.blastn.BlastnRunner.BLAST_TOOL;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.mapping.SequenceTester.FLD_SEQ_ID;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.BLASTN_WORD_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.MIN_BLAST_ALIGNMENT_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.blastn.BlastnMatch;
import com.hartwig.hmftools.common.blastn.BlastnRunner;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class BlastNSeqTester
{
    private final SeqTestConfig mConfig;
    private final BlastnRunner mBlastnRunner;

    public BlastNSeqTester(final SeqTestConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;

        if(configBuilder.hasValue(BLAST_TOOL) && configBuilder.hasValue(BLAST_DB))
        {
            mBlastnRunner = initialiseBlastNRunner(configBuilder.getValue(BLAST_TOOL), configBuilder.getValue(BLAST_DB));
        }
        else
        {
            mBlastnRunner = null;
        }
    }

    private BlastnRunner initialiseBlastNRunner(final String blastNTool, final String blastNDb)
    {
        if(blastNTool == null || blastNDb == null)
            return null;

        return new BlastnRunner.Builder()
                .withTask("megablast")
                .withPrefix(mConfig.OutputId)
                .withBlastDir(blastNTool)
                .withBlastDb(blastNDb)
                .withOutputDir(mConfig.OutputDir)
                .withKeepOutput(true)
                .withWordSize(BLASTN_WORD_SIZE)
                .withSubjectBestHit(true)
                .withNumThreads(mConfig.Threads)
                .build();
    }

    public void run(final List<SequenceInfo> sequenceInfos)
    {
        if(mBlastnRunner == null)
            return;

        long startTimeMs = System.currentTimeMillis();

        try
        {
            String outputFile = mConfig.OutputDir + mConfig.OutputId + ".blastn_results.tsv";

            BufferedWriter writer = initialiseBlastNWriter(outputFile);

            GU_LOGGER.info("running BlastN for {} sequences", sequenceInfos.size());

            Map<Integer,String> sequencesMap = Maps.newHashMap();

            for(SequenceInfo sequenceInfo : sequenceInfos)
            {
                sequencesMap.put(sequenceInfo.Id, sequenceInfo.Sequence);
            }

            Multimap<Integer, BlastnMatch> sequenceResults = mBlastnRunner.run(sequencesMap);

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

        GU_LOGGER.info("BlastN complete, mins({})", runTimeMinsStr(startTimeMs));
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

        SequenceInfo.addSequenceHeader(sj);

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

            // sequenceInfo.addSequenceInfo(sj);

            writer.write(sj.toString());
            writer.newLine();
            return;
        }

        List<BlastnMatch> sortedMatches = matches.stream()
                .filter(x -> isPrimaryBlastnMatch(x))
                .filter(x -> x.getAlignmentLength() >= MIN_BLAST_ALIGNMENT_LENGTH)
                .collect(Collectors.toList());

        Collections.sort(sortedMatches, Comparator.comparingDouble(x -> -x.BitScore));

        double sumBitScore = calcSumBitScore(sortedMatches, MIN_BLAST_ALIGNMENT_LENGTH);
        int resultCount = sortedMatches.size();

        // filter results as per usage for probes

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

            // sequenceInfo.addSequenceInfo(sj);

            writer.write(sj.toString());
            writer.newLine();
        }
    }

}
