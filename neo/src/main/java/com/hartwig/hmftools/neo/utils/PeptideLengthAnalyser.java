package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.stats.FisherExactTest;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class PeptideLengthAnalyser
{
    private final Map<String,List<BindPositionData>> mAlleleDataMap;
    private final Set<Integer> mPositions;
    private final Set<Integer> mPeptideLengths;
    private final String mOutputDir;
    private final BufferedWriter mWriter;

    private final FisherExactTest mFisherEt;

    private static final String NMER_FREQ_FILE = "nmer_freq_file";

    public PeptideLengthAnalyser(final CommandLine cmd)
    {
        String dataFile = cmd.getOptionValue(NMER_FREQ_FILE);

        mOutputDir = parseOutputDir(cmd);

        mWriter = initProbWriter();

        mAlleleDataMap = Maps.newHashMap();
        mPositions = Sets.newHashSet();
        mPeptideLengths = Sets.newHashSet();

        loadData(dataFile);

        mFisherEt = new FisherExactTest();
        mFisherEt.initialise(1000000);
    }

    public void run()
    {
        int maxPosition = mPositions.stream().mapToInt(x -> x.intValue()).max().orElse(0);

        for(Map.Entry<String,List<BindPositionData>> entry : mAlleleDataMap.entrySet())
        {
            String allele = entry.getKey();

            NE_LOGGER.info("evaluating allele({})", allele);

            for(int pos = 0; pos <= maxPosition; ++pos)
            {
                evaluate(allele, entry.getValue(), pos);
            }
        }

        closeBufferedWriter(mWriter);
    }

    private BufferedWriter initProbWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputDir + "peptide_length_probs.csv", false);

            writer.write("Allele,Position,AminoAcid,PeptideLength,Expected,Both,Prob");
            writer.write(",TotalBinds,WithLength,WithAminoAcid,WithAminoAcidNoLength,WithLengthNoAcid,Neither");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write study prob data: {}", e.toString());
        }

        return null;
    }

    private static final double LOW_PROB = 1e-6;

    private void evaluate(final String allele, final List<BindPositionData> posDataList, int position)
    {
        try
        {
            for(Character aminoAcid : AMINO_ACIDS)
            {
                for(Integer peptideLength : mPeptideLengths)
                {
                    int totalBinds = 0;
                    int withLength = 0;
                    int withAminoAcid = 0;
                    int withAminoAcidNoLength = 0;
                    int withLengthNoAcid = 0;
                    int both = 0;
                    int neither = 0;

                    for(BindPositionData data : posDataList)
                    {
                        if(data.Position != position)
                            continue;

                        totalBinds += data.BindCount;

                        if(data.AminoAcid == aminoAcid && data.PeptideLength == peptideLength)
                        {
                            withAminoAcid += data.BindCount;
                            withLength += data.BindCount;
                            both += data.BindCount;
                        }
                        else if(data.AminoAcid == aminoAcid && data.PeptideLength != peptideLength)
                        {
                            withAminoAcid += data.BindCount;
                            withAminoAcidNoLength += data.BindCount;
                        }
                        else if(data.AminoAcid != aminoAcid && data.PeptideLength == peptideLength)
                        {
                            withLengthNoAcid += data.BindCount;
                            withLength += data.BindCount;
                        }
                        else
                        {
                            neither += data.BindCount;
                        }
                    }

                    double expected = (withAminoAcid * withLength) / (double) totalBinds;

                    double prob = mFisherEt.calc(both, withAminoAcidNoLength, withLengthNoAcid, neither, expected);

                    if(prob < LOW_PROB)
                    {
                        NE_LOGGER.debug(String.format("allele(%s) position(%d) aminoAcids(%s) pepLen(%d) expected(%.3f) prob(%4.3e)",
                                allele, position, aminoAcid, peptideLength, expected, prob));

                        mWriter.write(String.format("%s,%s,%s,%d,%.2f,%d,%4.3e,%d,%d,%d,%d,%d,%d",
                                allele, position, aminoAcid, peptideLength, expected, both, prob, totalBinds,
                                withLength, withAminoAcid, withAminoAcidNoLength, withLengthNoAcid, neither));
                        mWriter.newLine();
                    }
                }
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write study prob data: {}", e.toString());
        }
    }

    private void loadData(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int pepLenIndex = fieldsIndexMap.get("PeptideLength");
            int posIndex = fieldsIndexMap.get("Position");
            int aaIndex = fieldsIndexMap.get("AminoAcid");
            int bcIndex = fieldsIndexMap.get("BindCount");

            String currentAllele = "";
            List<BindPositionData> currentBindList = null;

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                String allele = items[alleleIndex];

                if(!allele.equals(currentAllele))
                {
                    currentAllele = allele;
                    currentBindList = Lists.newArrayList();
                    mAlleleDataMap.put(allele, currentBindList);
                }

                BindPositionData data = new BindPositionData(
                        allele, Integer.parseInt(items[pepLenIndex]), items[aaIndex].charAt(0),
                        Integer.parseInt(items[posIndex]), Integer.parseInt(items[bcIndex]));

                mPositions.add(data.Position);
                mPeptideLengths.add(data.PeptideLength);

                currentBindList.add(data);
            }

            NE_LOGGER.info("loaded {} alleles, {} positions, {} peptide lengths with {} study data items from file({})",
                    mAlleleDataMap.size(), mPositions.size(), mPeptideLengths.size(),
                    mAlleleDataMap.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read study data file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(NMER_FREQ_FILE, true, "N-mer frequency data file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PeptideLengthAnalyser peptideLengthAnalyser = new PeptideLengthAnalyser(cmd);
        peptideLengthAnalyser.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private class BindPositionData
    {
        public final String Allele;
        public final int PeptideLength;
        public final Character AminoAcid;
        public final int Position;
        public final int BindCount;

        public BindPositionData(final String allele, final int peptideLength, final char aminoAcid, final int position, final int bindCount)
        {
            Allele = allele;
            PeptideLength = peptideLength;
            AminoAcid = aminoAcid;
            Position = position;
            BindCount = bindCount;
        }
    }
}
