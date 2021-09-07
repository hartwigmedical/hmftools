package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_DELIM;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.HG19;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.COMMON_ALLELES_FREQ_CUTOFF;
import static com.hartwig.hmftools.lilac.LilacConstants.EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_H;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_Y;
import static com.hartwig.hmftools.lilac.LilacConstants.STOP_LOSS_ON_C_ALLELE;
import static com.hartwig.hmftools.lilac.LilacConstants.getAminoAcidExonBoundaries;
import static com.hartwig.hmftools.lilac.LilacConstants.getNucleotideExonBoundaries;
import static com.hartwig.hmftools.lilac.hla.HlaContextFactory.populateNucleotideExonBoundaries;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.buildAminoAcidSequenceFromNucleotides;

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
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.lilac.cohort.CohortFrequency;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaAlleleCache;
import com.hartwig.hmftools.lilac.read.Indel;
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile;
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
    public final List<HlaSequenceLoci> HlaYAminoAcidSequences;

    public final List<HlaAllele> CommonAlleles; // common in population
    public final Map<Indel,HlaAllele> KnownStopLossIndelAlleles;

    private final CohortFrequency mAlleleFrequencies;

    public final Map<String,TranscriptData> HlaTranscriptData;

    public final LociPosition LociPositionFinder;

    private final HlaAlleleCache mAlleleCache;

    private HlaSequenceLoci mDeflatedSequenceTemplate;

    public static final String NUC_REF_FILE = "hla_ref_nucleotide_sequences.csv";
    public static final String AA_REF_FILE = "hla_ref_aminoacid_sequences.csv";

    private static final String COHORT_ALLELE_FREQ_FILE = "lilac_allele_frequencies.csv";

    // sequence used to printing amino acid sequences to file
    public static final HlaAllele DEFLATE_TEMPLATE = HlaAllele.fromString("A*01:01");

    public static Indel STOP_LOSS_ON_C_INDEL = null;

    public static final List<Indel> INDEL_PON = Lists.newArrayList();

    public ReferenceData(final String resourceDir, final LilacConfig config)
    {
        mResourceDir = resourceDir;
        mConfig = config;

        mAlleleCache = new HlaAlleleCache();

        mAlleleFrequencies = new CohortFrequency(!mResourceDir.isEmpty() ? mResourceDir + COHORT_ALLELE_FREQ_FILE : "");

        NucleotideSequences = Lists.newArrayList();
        AminoAcidSequences = Lists.newArrayList();
        AminoAcidSequencesWithInserts = Lists.newArrayList();
        AminoAcidSequencesWithDeletes = Lists.newArrayList();
        HlaYNucleotideSequences = Lists.newArrayList();
        HlaYAminoAcidSequences = Lists.newArrayList();

        mDeflatedSequenceTemplate = null;

        setKnownStopLossIndels(config.RefGenVersion);
        setPonIndels(config.RefGenVersion);

        HlaTranscriptData = Maps.newHashMap();

        // load gene definitions and other constants
        populateHlaTranscripts(HlaTranscriptData, config.RefGenVersion);
        LociPositionFinder = new LociPosition(HlaTranscriptData.values().stream().collect(Collectors.toList()));

        populateNucleotideExonBoundaries();

        CommonAlleles = Lists.newArrayList();
        KnownStopLossIndelAlleles = Maps.newHashMap();
    }

    private static void setPonIndels(final RefGenomeVersion version)
    {
        // load indel PON
        String refFile = version.is37() ? "/pon/indels_v37.csv" : "/pon/indels_v38.csv";

        final List<String> ponLines = new BufferedReader(new InputStreamReader(
                RefGenomeCoordinates.class.getResourceAsStream(refFile)))
                .lines().collect(Collectors.toList());

        ponLines.stream().map(x -> Indel.fromString(x)).forEach(x -> INDEL_PON.add(x));
    }

    private static void setKnownStopLossIndels(final RefGenomeVersion version)
    {
        int position = version.is37() ? 31237115 : 31269338;
        STOP_LOSS_ON_C_INDEL = new Indel("6", position, "CN", "C");

        if(version == V38 || version == HG19)
            LilacConstants.HLA_CHR = "chr6"; // should be abe to get from transcript info at time of use when get rid of NamedBed
    }

    public CohortFrequency getAlleleFrequencies() { return mAlleleFrequencies; }

    public HlaAllele findAllele(final String alleleStr, boolean isFourDigit)
    {
        return isFourDigit ? mAlleleCache.findFourDigitAllele(alleleStr) : mAlleleCache.findAllele(alleleStr);
    }

    public boolean load()
    {
        if(!mResourceDir.isEmpty())
        {
            String nucleotideFilename = mResourceDir + NUC_REF_FILE;

            LL_LOGGER.info("reading nucleotide file: {}", nucleotideFilename);

            if(!loadSequenceFile(nucleotideFilename, NucleotideSequences, false))
                return false;

            String aminoAcidFilename = mResourceDir + AA_REF_FILE;

            LL_LOGGER.info("reading protein file: {}", aminoAcidFilename);

            if(!loadSequenceFile(aminoAcidFilename, AminoAcidSequences, true))
                return false;
        }

        // load and register configured and known alleles
        mAlleleCache.rebuildProteinAlleles(mConfig.ActualAlleles);
        mAlleleCache.rebuildProteinAlleles(mConfig.RestrictedAlleles);

        // apply PON
        // "A*01:81", "A*01:237", "A*11:126", "A*11:353", "A*25:68", "A*30:95", "A*30:136", "A*31:135", "A*33:191");

        loadCommonAlleles();

        if(!CommonAlleles.isEmpty())
        {
            LL_LOGGER.info("loaded {} common alleles", CommonAlleles.size());
        }

        loadStopLossRecoveryAllele();

        HlaYNucleotideSequences.addAll(NucleotideSequences.stream().filter(x -> x.Allele.Gene.equals(GENE_Y)).collect(Collectors.toList()));
        HlaYNucleotideSequences.forEach(x -> NucleotideSequences.remove(x));

        for(HlaSequenceLoci sequenceLoci : AminoAcidSequences)
        {
            if(sequenceLoci.hasInserts())
            {
                AminoAcidSequencesWithInserts.add(sequenceLoci);
            }
            else if(sequenceLoci.hasDeletes() || KnownStopLossIndelAlleles.values().stream().anyMatch(x -> x.equals(sequenceLoci.Allele)))
            {
                AminoAcidSequencesWithDeletes.add(sequenceLoci);
            }
        }

        markExonBoundaryWildcards();
        buildHlaYAminoAcidSequences();

        return true;
    }

    private void buildHlaYAminoAcidSequences()
    {
        HlaYAminoAcidSequences.addAll(AminoAcidSequences.stream().filter(x -> x.Allele.Gene.equals(GENE_Y)).collect(Collectors.toList()));

        if(!HlaYAminoAcidSequences.isEmpty())
        {
            HlaYAminoAcidSequences.forEach(x -> AminoAcidSequences.remove(x));
        }
        else
        {
            for(HlaSequenceLoci sequenceLoci : HlaYNucleotideSequences)
            {
                HlaYAminoAcidSequences.add(buildAminoAcidSequenceFromNucleotides(sequenceLoci, mDeflatedSequenceTemplate));
            }
        }
    }

    private void markExonBoundaryWildcards()
    {
        for(HlaSequenceLoci sequence : AminoAcidSequences)
        {
            List<Integer> exonBoundaries = getAminoAcidExonBoundaries(sequence.Allele.Gene);
            sequence.setExonBoundaryWildcardsWildcards(exonBoundaries);
        }

        for(HlaSequenceLoci sequence : NucleotideSequences)
        {
            List<Integer> exonBoundaries = getNucleotideExonBoundaries(sequence.Allele.Gene);
            sequence.setExonBoundaryWildcardsWildcards(exonBoundaries);
        }
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
        KnownStopLossIndelAlleles.put(STOP_LOSS_ON_C_INDEL, mAlleleCache.requestFourDigit(STOP_LOSS_ON_C_ALLELE));
    }

    private boolean excludeAllele(final HlaAllele allele)
    {
        if(allele.Gene.equals(GENE_H))
            return true;

        final HlaAllele allele4d = allele.asFourDigit();

        if(EXCLUDED_ALLELES.stream().anyMatch(x -> allele4d.matches(x)))
            return true;

        if(mConfig == null)
            return false;

        if(!mConfig.RestrictedAlleles.isEmpty())
        {
            if(mConfig.ActualAlleles.stream().anyMatch(x -> x.matches(allele4d)))
                return false;

            if(mConfig.RestrictedAlleles.stream().noneMatch(x -> x.matches(allele4d)))
                return true;
        }

        return false;
    }

    public static void populateHlaTranscripts(final Map<String,TranscriptData> hlaTranscriptMap, final RefGenomeVersion version)
    {
        String transcriptsFile = version.is37() ? "/alleles/hla_transcripts_v37.csv" : "/alleles/hla_transcripts_v38.csv";

        final List<String> hlaTranscriptData = new BufferedReader(new InputStreamReader(
                ReferenceData.class.getResourceAsStream(transcriptsFile))).lines().collect(Collectors.toList());

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
            fileContents.remove(0); // remove header
            loadSequenceFile(fileContents, sequenceData, isProteinFile);
            LL_LOGGER.info("loaded {} sequences from file {}", sequenceData.size(), filename);
            return true;
        }
        catch (IOException e)
        {
            LL_LOGGER.error("failed to load ref sequence data from file({}): {}", filename, e.toString());
            return false;
        }
    }

    public void loadSequenceFile(final List<String> fileContents, final List<HlaSequenceLoci> sequenceData, boolean isProteinFile)
    {
        for(String line : fileContents)
        {
            String[] items = line.split(",");

            if(items.length != 2)
                return;

            String alleleStr = items[0];

            HlaAllele allele = isProteinFile ? mAlleleCache.requestFourDigit(alleleStr) : mAlleleCache.request(alleleStr);

            boolean isDefaultTemplate = isProteinFile && mDeflatedSequenceTemplate == null && allele.matches(DEFLATE_TEMPLATE);
            boolean excludeAllele = excludeAllele(allele);

            if(!isDefaultTemplate && excludeAllele)
                continue;

            String sequenceStr = items[1];

            HlaSequenceLoci newSequence = HlaSequenceFile.createFromReference(allele, sequenceStr, isProteinFile);

            if(isDefaultTemplate)
                mDeflatedSequenceTemplate = newSequence;

            if(!excludeAllele)
                sequenceData.add(newSequence);
        }
    }
}

