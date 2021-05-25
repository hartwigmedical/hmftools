package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.COMMON_ALLELES_FREQ_CUTOFF;
import static com.hartwig.hmftools.lilac.LilacConstants.EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_Y;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.lilac.cohort.CohortFrequency;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaAlleleCache;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class ReferenceData
{
    private final String mResourceDir;
    private final LilacConfig mConfig;

    public final List<HlaSequenceLoci> NucleotideSequences;
    public final List<HlaSequenceLoci> AminoAcidSequences;
    public final List<HlaSequenceLoci> AminoAcidSequencesWithInserts;
    public final List<HlaSequenceLoci> AminoAcidSequencesWithDeletes;
    public final List<HlaSequenceLoci> HlaYNucleotideSequences;

    public final List<HlaAllele> CommonAlleles; // common in population
    public final List<HlaAllele> StopLossRecoveryAlleles;

    private final CohortFrequency mAlleleFrequencies;

    public final Map<String,TranscriptData> HlaTranscriptData;

    public final LociPosition LociPositionFinder;

    private final HlaAlleleCache mAlleleCache;

    private HlaSequenceLoci mDeflatedSequenceTemplate;

    public static final char SEQUENCE_DELIM = '|';
    public static final String NUC_REF_FILE = "hla_ref_nucleotide_sequences.csv";
    public static final String AA_REF_FILE = "hla_ref_aminoacid_sequences.csv";

    private static final String COHORT_ALLELE_FREQ_FILE = "lilac_allele_frequencies.csv";

    // sequence used to printing amino acid sequences to file
    public static final HlaAllele DEFLATE_TEMPLATE = HlaAllele.fromString("A*01:01");

    public ReferenceData(final String resourceDir, final LilacConfig config)
    {
        mResourceDir = resourceDir;
        mConfig = config;

        mAlleleCache = new HlaAlleleCache();

        mAlleleFrequencies = new CohortFrequency(mResourceDir + COHORT_ALLELE_FREQ_FILE);

        NucleotideSequences = Lists.newArrayList();
        AminoAcidSequences = Lists.newArrayList();
        AminoAcidSequencesWithInserts = Lists.newArrayList();
        AminoAcidSequencesWithDeletes = Lists.newArrayList();
        HlaYNucleotideSequences = Lists.newArrayList();

        mDeflatedSequenceTemplate = null;

        HlaTranscriptData = Maps.newHashMap();

        // load gene definitions and other constants
        populateHlaTranscripts(HlaTranscriptData);
        LociPositionFinder = new LociPosition(HlaTranscriptData.values().stream().collect(Collectors.toList()));

        CommonAlleles = Lists.newArrayList();
        StopLossRecoveryAlleles = Lists.newArrayList();
    }

    public HlaSequenceLoci getDeflatedSequenceTemplate() { return mDeflatedSequenceTemplate; }
    public CohortFrequency getAlleleFrequencies() { return mAlleleFrequencies; }

    public boolean load()
    {
        // load and register configured and known alleles
        mAlleleCache.rebuildProteinAlleles(mConfig.ExpectedAlleles);
        mAlleleCache.rebuildProteinAlleles(mConfig.RestrictedAlleles);

        loadCommonAlleles();

        if(!CommonAlleles.isEmpty())
        {
            LL_LOGGER.info("loaded {} common alleles", CommonAlleles.size());
        }

        loadStopLossRecoveryAllele();

        LL_LOGGER.info("reading nucleotide files");

        String nucleotideFilename = mResourceDir + NUC_REF_FILE;

        if(!loadSequenceFile(nucleotideFilename, NucleotideSequences, false))
            return false;

        HlaYNucleotideSequences.addAll(NucleotideSequences.stream().filter(x -> x.Allele.Gene.equals(GENE_Y)).collect(Collectors.toList()));
        HlaYNucleotideSequences.forEach(x -> NucleotideSequences.remove(x));

        LL_LOGGER.info("reading protein files");

        String aminoAcidFilename = mResourceDir + AA_REF_FILE;

        if(!loadSequenceFile(aminoAcidFilename, AminoAcidSequences, true))
            return false;

        AminoAcidSequencesWithInserts.addAll(AminoAcidSequences.stream().filter(x -> x.containsInserts()).collect(Collectors.toList()));
        AminoAcidSequencesWithDeletes.addAll(AminoAcidSequences.stream().filter(x -> x.containsDeletes()).collect(Collectors.toList()));

        return true;
    }

    private void loadCommonAlleles()
    {
        mAlleleFrequencies.getAlleleFrequencies().entrySet().stream()
                .filter(x -> x.getValue() >= COMMON_ALLELES_FREQ_CUTOFF)
                .map(x -> mAlleleCache.requestFourDigit(x.getKey().toString()))
                .forEach(x -> CommonAlleles.add(x));
    }

    private void loadStopLossRecoveryAllele()
    {
        StopLossRecoveryAlleles.add(mAlleleCache.requestFourDigit("C*04:09N"));
    }

    public List<HlaAllele> getWildcardAlleles()
    {
        return AminoAcidSequences.stream().filter(x -> x.containsWildcards()).map(x -> x.Allele).collect(Collectors.toList());
    }

    private boolean excludeAllele(final HlaAllele allele)
    {
        final HlaAllele allele4d = allele.asFourDigit();

        if(EXCLUDED_ALLELES.stream().anyMatch(x -> allele4d.matches(x)))
            return true;

        if(mConfig == null)
            return false;

        if(!mConfig.RestrictedAlleles.isEmpty())
        {
            if(mConfig.ExpectedAlleles.stream().anyMatch(x -> x.matches(allele4d)))
                return false;

            if(mConfig.RestrictedAlleles.stream().noneMatch(x -> x.matches(allele4d)))
                return true;
        }

        return false;
    }

    public static void populateHlaTranscripts(final Map<String,TranscriptData> hlaTranscriptMap)
    {
        final List<String> hlaTranscriptData = new BufferedReader(new InputStreamReader(
                ReferenceData.class.getResourceAsStream("/alleles/hla_transcripts_v37.csv")))
                .lines().collect(Collectors.toList());

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(hlaTranscriptData.get(0), ENSEMBL_DELIM);
        hlaTranscriptData.remove(0);

        int geneIdIndex = fieldsIndexMap.get("GeneId");
        int geneNameIndex = fieldsIndexMap.get("GeneName");
        int strandIndex = fieldsIndexMap.get("Strand");
        int transIdIndex = fieldsIndexMap.get("TransId");
        int transNameIndex = fieldsIndexMap.containsKey("TransName") ? fieldsIndexMap.get("TransName") : fieldsIndexMap.get("Trans");
        int biotypeIndex = fieldsIndexMap.get("BioType");
        int transStartIndex = fieldsIndexMap.get("TransStart");
        int transEndIndex = fieldsIndexMap.get("TransEnd");
        int exonRankIndex = fieldsIndexMap.get("ExonRank");
        int exonStartIndex = fieldsIndexMap.get("ExonStart");
        int exonEndIndex = fieldsIndexMap.get("ExonEnd");
        int exonPhaseIndex = fieldsIndexMap.get("ExonPhase");
        int exonEndPhaseIndex = fieldsIndexMap.get("ExonEndPhase");
        int codingStartIndex = fieldsIndexMap.get("CodingStart");
        int codingEndIndex = fieldsIndexMap.get("CodingEnd");

        String currentGene = "";
        TranscriptData currentTrans = null;
        List<ExonData> exonDataList = null;

        for(String line : hlaTranscriptData)
        {
            String[] items = line.split(ENSEMBL_DELIM);

            String geneId = items[geneIdIndex];
            int transId = Integer.parseInt(items[transIdIndex]);

            if(!geneId.equals(currentGene))
            {
                currentGene = geneId;

                String geneName = items[geneNameIndex];

                Integer codingStart = Integer.parseInt(items[codingStartIndex]);
                Integer codingEnd = Integer.parseInt(items[codingEndIndex]);

                currentTrans = new TranscriptData(
                        transId, items[transNameIndex], geneId, true, Byte.parseByte(items[strandIndex]),
                        Integer.parseInt(items[transStartIndex]), Integer.parseInt(items[transEndIndex]),
                        codingStart, codingEnd, items[biotypeIndex]);

                hlaTranscriptMap.put(geneName, currentTrans);

                exonDataList = currentTrans.exons();
            }

            ExonData exonData = new ExonData(
                    transId, Integer.parseInt(items[exonStartIndex]), Integer.parseInt(items[exonEndIndex]),
                    Integer.parseInt(items[exonRankIndex]), Integer.parseInt(items[exonPhaseIndex]), Integer.parseInt(items[exonEndPhaseIndex]));

            exonDataList.add(exonData);
        }
    }

    private boolean loadSequenceFile(final String filename, final List<HlaSequenceLoci> sequenceData, boolean isProteinFile)
    {
        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            fileContents.remove(0);

            for(String line : fileContents)
            {
                String[] items = line.split(",");

                if(items.length != 2)
                    return false;

                String alleleStr = items[0];
                HlaAllele allele = isProteinFile ? mAlleleCache.requestFourDigit(alleleStr) : mAlleleCache.request(alleleStr);

                boolean isDefaultTemplate = isProteinFile && mDeflatedSequenceTemplate == null && allele.matches(DEFLATE_TEMPLATE);
                boolean excludeAllele = excludeAllele(allele);

                if(!isDefaultTemplate && excludeAllele)
                    continue;

                List<String> sequences = Lists.newArrayList();
                String sequenceStr = items[1];

                int index = 0;
                String sequence = "";
                boolean inMulti = false;
                while(index < sequenceStr.length())
                {
                    char nextChar = sequenceStr.charAt(index);
                    boolean isMulti = nextChar == SEQUENCE_DELIM;

                    if(inMulti || isMulti)
                    {
                        if(inMulti && isMulti)
                        {
                            inMulti = false;
                            sequences.add(sequence);
                            sequence = "";
                        }
                        else if(isMulti)
                        {
                            // start of new multi-char sequence
                            inMulti = true;
                        }
                        else
                        {
                            sequence += nextChar;
                        }
                    }
                    else
                    {
                        sequences.add(String.valueOf(nextChar));
                    }

                    ++index;
                }

                HlaSequenceLoci newSequence = new HlaSequenceLoci(allele, sequences);

                if(isDefaultTemplate)
                    mDeflatedSequenceTemplate = newSequence;

                if(!excludeAllele)
                    sequenceData.add(newSequence);
            }

            LL_LOGGER.info("loaded {} sequences from file {}", sequenceData.size(), filename);
            return true;
        }
        catch (IOException e)
        {
            LL_LOGGER.error("failed to load ref sequence data from file({}): {}", filename, e.toString());
            return false;
        }
    }
}

