package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELES;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PATIENT_ID;
import static com.hartwig.hmftools.neo.bind.BindCommon.ITEM_DELIM;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTHS;
import static com.hartwig.hmftools.neo.bind.BinderConfig.REQUIRED_PEPTIDE_LENGTHS;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_BASE_COUNT;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.utils.FileWriterUtils;
import com.hartwig.hmftools.neo.bind.RandomPeptideConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class GartnerDataPrep
{
    private final String mOutputDir;
    private final List<Integer> mRequiredPeptideLengths;

    private final Map<String,List<TranscriptAminoAcids>> mGeneAminoAcidMap;

    private final Map<String,List<String>> mPatientAlleles;
    private final List<MutationData> mMutations;

    // config
    private static final String PATIENT_ALLELES_FILE = "patient_alleles_file";
    private static final String MUTATIONS_FILE = "mutations_file";

    public GartnerDataPrep(final CommandLine cmd)
    {
        mPatientAlleles = Maps.newHashMap();
        mMutations = Lists.newArrayList();

        String ensemblDataDir = cmd.getOptionValue(ENSEMBL_DATA_DIR);
        mGeneAminoAcidMap = Maps.newHashMap();

        Map<String,TranscriptAminoAcids> transAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, transAminoAcidMap, Lists.newArrayList());

        for(TranscriptAminoAcids transAminoAcids : transAminoAcidMap.values())
        {
            List<TranscriptAminoAcids> genesList = mGeneAminoAcidMap.get(transAminoAcids.GeneName);
            if(genesList == null)
            {
                genesList = Lists.newArrayList();
                mGeneAminoAcidMap.put(transAminoAcids.GeneName, genesList);
            }

            genesList.add(transAminoAcids);
        }

        loadPatientAlleles(cmd.getOptionValue(PATIENT_ALLELES_FILE));
        loadMutations(cmd.getOptionValue(MUTATIONS_FILE));

        mRequiredPeptideLengths = DEFAULT_PEPTIDE_LENGTHS;

        mOutputDir = FileWriterUtils.parseOutputDir(cmd);
    }

    public void run()
    {
        if(mPatientAlleles.isEmpty() || mMutations.isEmpty())
        {
            System.exit(1);
            return;
        }

        BufferedWriter writer = initWriter();

        for(MutationData mutation : mMutations)
        {
            generatePeptides(mutation);
            writeMutation(writer, mutation);
        }

        closeBufferedWriter(writer);
    }

    private boolean generatePeptides(final MutationData mutation)
    {
        String mutationAAs = mutation.MutEpitope;
        String wildtypeAAs = mutation.WtEpitope;

        int diffIndex = -1;

        for(int i = 0; i < wildtypeAAs.length(); ++i)
        {
            if(i >= mutationAAs.length())
                break;

            if(wildtypeAAs.charAt(i) != mutationAAs.charAt(i))
            {
                diffIndex = i;
                break;
            }
        }

        if(diffIndex < 0)
        {
            NE_LOGGER.warn("mutation({}) no diff in AAs", mutation);
            return false;
        }

        boolean sameAaCount = mutationAAs.length() == wildtypeAAs.length();

        for(int i = 0; i < mRequiredPeptideLengths.size(); ++i)
        {
            int peptideLength = mRequiredPeptideLengths.get(i);

            int startIndex = diffIndex - peptideLength + 1;

            if(startIndex < 0)
                continue;

            int endDiffIndex = sameAaCount ? diffIndex : mutationAAs.length() - 1;

            for(; startIndex <= endDiffIndex; ++startIndex)
            {
                int endIndex = startIndex + peptideLength;
                if(endIndex >= mutationAAs.length())
                    break;

                String peptide = mutationAAs.substring(startIndex, endIndex);

                if(wildtypeAAs.contains(peptide)) // double-check doesn't revert to wildtype bases,
                    continue;

                String upFlank = "";
                String downFlank = "";

                int upFlankIndex = max(startIndex - FLANK_BASE_COUNT, 0);
                if(upFlankIndex < startIndex)
                    upFlank = mutationAAs.substring(upFlankIndex, startIndex);

                int downFlankIndex = min(endIndex + FLANK_BASE_COUNT, mutationAAs.length() - 1);
                if(downFlankIndex > endIndex)
                    downFlank = mutationAAs.substring(endIndex, downFlankIndex);

                // String[] flanks = getFlanks(mutation.GeneName, peptide);

                mutation.Peptides.add(new PeptideData(peptide, upFlank, downFlank));
            }
        }

        return true;
    }

    private String[] getFlanks(final String geneName, final String peptide)
    {
        String[] flanks = new String[] { "", ""};

        List<TranscriptAminoAcids> genesList = mGeneAminoAcidMap.get(geneName);
        if(genesList == null)
            return flanks;

        for(TranscriptAminoAcids transAminoAcids : genesList)
        {
            int aaIndex = transAminoAcids.AminoAcids.indexOf(peptide);
            if(aaIndex < 0)
                continue;

            int upFlankBases = min(aaIndex, FLANK_BASE_COUNT);

            if(upFlankBases > 0)
            {
                flanks[0] = transAminoAcids.AminoAcids.substring(aaIndex - upFlankBases, aaIndex);
            }

            int downFlankStartPos = aaIndex + peptide.length();
            int downFlankBases = min(transAminoAcids.AminoAcids.length() - downFlankStartPos - 1, FLANK_BASE_COUNT);

            if(downFlankBases > 0)
            {
                flanks[1] = transAminoAcids.AminoAcids.substring(downFlankStartPos, downFlankStartPos + downFlankBases);
            }

            break;
        }

        return flanks;
    }

    private BufferedWriter initWriter()
    {
        try
        {
            String filename = mOutputDir + "gart_allele_peptide_data";

            BufferedWriter writer = createBufferedWriter(filename, false);
            writer.write("PatientId,Allele,Peptide,UpFlank,DownFlank,OtherInfo");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write allele peptide data: {}", e.toString());
            return null;
        }
    }

    private void writeMutation(final BufferedWriter writer, final MutationData mutation)
    {
        List<String> alleles = mPatientAlleles.get(mutation.PatientId);

        if(alleles == null)
        {
            NE_LOGGER.error("mutation({}) missing allele data by patientId", mutation);
            return;
        }

        try
        {
            for(String allele : alleles)
            {
                for(PeptideData peptideData : mutation.Peptides)
                {
                    writer.write(String.format("%s,%s,%s,%s,%s,%s_%s",
                            mutation.PatientId, allele, peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank,
                            mutation.GeneName, mutation.VariantKey));

                    writer.newLine();
                }
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write mutation allele peptide data: {}", e.toString());
        }
    }

    private boolean loadPatientAlleles(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            // PatientId,Alleles
            String header = lines.get(0);
            lines.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELES);
            int patientIndex = fieldsIndexMap.get(FLD_PATIENT_ID);

            for(String line : lines)
            {
                String[] values = line.split(DELIM, -1);
                String patient = values[patientIndex];
                String allelesStr = values[alleleIndex];

                List<String> alleles = Arrays.stream(allelesStr.split(ITEM_DELIM)).collect(Collectors.toList());

                mPatientAlleles.put(patient, alleles);
            }

            NE_LOGGER.info("loaded {} patient alleles", mPatientAlleles.size());
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load patient-alleles file({}): {}", filename, e.toString());
        }

        return true;
    }

    private boolean loadMutations(final String filename)
    {
        // VariantKey,PatientId,GeneName,MutationType,AaChange,WtEpitope,MutEpitope
        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            // PatientId,Alleles
            String header = lines.get(0);
            lines.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);

            for(String line : lines)
            {
                mMutations.add(new MutationData(fieldsIndexMap, line));
            }

            NE_LOGGER.info("loaded {} mutations", mMutations.size());
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load mutations from file({}): {}", filename, e.toString());
        }

        return true;
    }

    private class PeptideData
    {
        public final String Peptide;
        public final String UpFlank;
        public final String DownFlank;

        public PeptideData(final String peptide, final String upFlank, final String downFlank)
        {
            Peptide = peptide;
            UpFlank = upFlank;
            DownFlank = downFlank;
        }
    }

    private class MutationData
    {
        // VariantKey,PatientId,GeneName,MutationType,AaChange,WtEpitope,MutEpitope
        public final String VariantKey;
        public final String PatientId;
        public final String GeneName;
        public final String MutationType;
        public final String AaChange;
        public final String WtEpitope;
        public final String MutEpitope;

        public final List<PeptideData> Peptides;

        public MutationData(final Map<String,Integer> fieldsIndexMap, final String line)
        {
            String[] values = line.split(DELIM, -1);

            VariantKey = values[fieldsIndexMap.get("VariantKey")];
            PatientId = values[fieldsIndexMap.get("PatientId")];
            GeneName = values[fieldsIndexMap.get("GeneName")];
            MutationType = values[fieldsIndexMap.get("MutationType")];
            AaChange = values[fieldsIndexMap.get("AaChange")];
            WtEpitope = values[fieldsIndexMap.get("WtEpitope")];
            MutEpitope = values[fieldsIndexMap.get("MutEpitope")];

            Peptides = Lists.newArrayList();
        }

        public String toString() { return String.format("patient(%s) var(%s) gene(%s)", PatientId, VariantKey, GeneName); }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        RandomPeptideConfig.addCmdLineArgs(options);
        options.addOption(PATIENT_ALLELES_FILE, true, "MCF predictions file");
        options.addOption(MUTATIONS_FILE, true, "Binding validation file");
        options.addOption(REQUIRED_PEPTIDE_LENGTHS, true, "Peptide lengths");
        addEnsemblDir(options);
        addLoggingOptions(options);
        addOutputDir(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        GartnerDataPrep gartnerDataPrep = new GartnerDataPrep(cmd);
        gartnerDataPrep.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
