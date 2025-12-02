package com.hartwig.hmftools.lilac;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_DELIM;
import static com.hartwig.hmftools.common.hla.HlaCommon.HLA_CHROMOSOME_V37;
import static com.hartwig.hmftools.common.hla.HlaCommon.HLA_CHROMOSOME_V38;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.CLASS_1_EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.LilacConstants.COMMON_ALLELES_FREQ_CUTOFF;
import static com.hartwig.hmftools.lilac.LilacConstants.CURRENT_GENES;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_CHR;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_DRB1_EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.LilacConstants.STOP_LOSS_ON_C_ALLELE;
import static com.hartwig.hmftools.lilac.LilacConstants.V37_HLA_REGION;
import static com.hartwig.hmftools.lilac.LilacConstants.V38_HLA_REGION;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_H;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_Y;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.buildAminoAcidSequenceFromNucleotides;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaAlleleCache;
import com.hartwig.hmftools.lilac.hla.HlaContextFactory;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.Indel;
import com.hartwig.hmftools.lilac.seq.HlaExonSequences;
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ReferenceData
{
    private final String mResourceDir;
    private final LilacConfig mConfig;

    private static String HLA_REF_BASES = null;
    private static int HLA_REF_START = 0;

    public static GeneCache GENE_CACHE = null;

    public final List<HlaSequenceLoci> NucleotideSequences;
    public final List<HlaSequenceLoci> AminoAcidSequences;
    public final List<HlaSequenceLoci> AminoAcidSequencesWithInserts;
    public final List<HlaSequenceLoci> AminoAcidSequencesWithDeletes;
    public final List<HlaSequenceLoci> HlaYNucleotideSequences;
    public final List<HlaSequenceLoci> HlaYAminoAcidSequences;

    // four-digit allele to seq
    public final Map<HlaAllele, HlaExonSequences> ExonSequencesLookup;

    public final List<HlaAllele> CommonAlleles; // common in population
    public final Map<Indel, HlaAllele> KnownStopLossIndelAlleles;

    private final CohortFrequency mAlleleFrequencies;

    // temporary until HlaContextFactor and are refactored to be unaware of class-type
    public static Map<HlaGene, List<Integer>> GENE_EXON_BOUNDARIES = null;

    public static HlaContextFactory HLA_CONTEXT_FACTORY = null;
    public static NucleotideGeneEnrichment NUC_GENE_FRAG_ENRICHMENT = null;

    private final HlaAlleleCache mAlleleCache;

    private HlaSequenceLoci mDeflatedSequenceTemplate;

    // external reference files
    public static final String NUC_REF_FILE = "hla_ref_nucleotide_sequences.csv";
    public static final String AA_REF_FILE = "hla_ref_aminoacid_sequences.csv";
    private static final String COHORT_ALLELE_FREQ_FILE = "lilac_allele_frequencies.csv";

    // sequence used to printing amino acid sequences to file
    public static final HlaAllele DEFLATE_TEMPLATE = HlaAllele.fromString("A*01:01");

    private static final List<String> EXCLUDED_ALLELES = Lists.newArrayList();

    public static Indel STOP_LOSS_ON_C_INDEL = null;

    public static List<Indel> INDEL_PON = null;

    public ReferenceData(final String resourceDir, final LilacConfig config)
    {
        mResourceDir = resourceDir;
        mConfig = config;

        Map<HlaGene, TranscriptData> hlaTranscriptMap = loadHlaTranscripts(config.RefGenVersion, CURRENT_GENES);
        if(CURRENT_GENES == GeneSelector.HLA_DRB1)
            trimDrb1Transcripts(hlaTranscriptMap);

        HLA_CHR = config.RefGenVersion.is38() ? HLA_CHROMOSOME_V38 : HLA_CHROMOSOME_V37;

        if(!"".equals(mConfig.RefGenome) && HLA_REF_BASES == null)
        {
            IndexedFastaSequenceFile refFastaSeqFile = null;
            try
            {
                refFastaSeqFile = new IndexedFastaSequenceFile(new File(mConfig.RefGenome));
            }
            catch(FileNotFoundException e)
            {
                LL_LOGGER.error("Failed to load ref genome: {}", mConfig.RefGenome);
                System.exit(1);
            }

            RefGenomeSource refGenome = new RefGenomeSource(refFastaSeqFile);
            if(config.RefGenVersion.is38())
            {
                HLA_REF_START = V38_HLA_REGION.start();
                HLA_REF_BASES = refGenome.getBaseString(V38_HLA_REGION.Chromosome, V38_HLA_REGION.start(), V38_HLA_REGION.end());
            }
            else
            {
                HLA_REF_START = V37_HLA_REGION.start();
                HLA_REF_BASES = refGenome.getBaseString(V37_HLA_REGION.Chromosome, V37_HLA_REGION.start(), V37_HLA_REGION.end());
            }
        }

        GENE_CACHE = new GeneCache(hlaTranscriptMap);
        EXCLUDED_ALLELES.clear();
        if(CURRENT_GENES.coversMhcClass1())
            EXCLUDED_ALLELES.addAll(CLASS_1_EXCLUDED_ALLELES);

        if(CURRENT_GENES.contains(HlaGene.HLA_DRB1))
            EXCLUDED_ALLELES.addAll(HLA_DRB1_EXCLUDED_ALLELES);

        // see note above
        GENE_EXON_BOUNDARIES = Maps.newHashMap();
        for(HlaGene gene : GENE_CACHE.GeneNames)
            GENE_EXON_BOUNDARIES.computeIfAbsent(gene, k -> Lists.newArrayList()).addAll(GENE_CACHE.AminoAcidExonBoundaries.get(gene));

        HLA_CONTEXT_FACTORY = new HlaContextFactory(GENE_EXON_BOUNDARIES);
        NUC_GENE_FRAG_ENRICHMENT = NucleotideGeneEnrichment.create(GENE_EXON_BOUNDARIES);

        mAlleleCache = new HlaAlleleCache();

        mAlleleFrequencies = new CohortFrequency(CURRENT_GENES, !mResourceDir.isEmpty() ? mResourceDir + COHORT_ALLELE_FREQ_FILE : "");

        NucleotideSequences = Lists.newArrayList();
        AminoAcidSequences = Lists.newArrayList();
        AminoAcidSequencesWithInserts = Lists.newArrayList();
        AminoAcidSequencesWithDeletes = Lists.newArrayList();
        HlaYNucleotideSequences = CURRENT_GENES.coversMhcClass1() ? Lists.newArrayList() : null;
        HlaYAminoAcidSequences = CURRENT_GENES.coversMhcClass1() ? Lists.newArrayList() : null;

        ExonSequencesLookup = Maps.newHashMap();

        mDeflatedSequenceTemplate = null;

        setPonIndels(config.RefGenVersion);

        CommonAlleles = Lists.newArrayList();
        KnownStopLossIndelAlleles = Maps.newHashMap();
    }

    public static String refBases(int startPos, int endPos)
    {
        return HLA_REF_BASES.substring(startPos - HLA_REF_START, endPos - HLA_REF_START + 1);
    }

    private static void trimDrb1Transcripts(final Map<HlaGene, TranscriptData> transcripts)
    {
        // restrict to first three exons
        TranscriptData transcript = transcripts.get(HlaGene.HLA_DRB1);
        List<ExonData> exons = Lists.newArrayList(transcript.exons());
        exons = exons.subList(exons.size() - 3, exons.size());

        // adjust end of exon 3 so that we get a multiple of 3 nucleotide coding bases
        int nucCount = 0;
        for(ExonData exon : exons)
        {
            int start = max(transcript.CodingStart, exon.Start);
            int end = min(transcript.CodingEnd, exon.End);
            nucCount += end - start + 1;
        }

        int trimmedCodingStart = exons.get(0).Start + (nucCount % 3);
        TranscriptData trimmedTranscript =
                new TranscriptData(transcript.TransId, transcript.TransName, transcript.GeneId, transcript.IsCanonical, transcript.Strand, transcript.TransStart, transcript.TransEnd, trimmedCodingStart, transcript.CodingEnd, transcript.BioType, transcript.RefSeqId);
        trimmedTranscript.setExons(exons);
        transcripts.put(HlaGene.HLA_DRB1, trimmedTranscript);
    }

    private static void setPonIndels(final RefGenomeVersion version)
    {
        if(INDEL_PON != null)
            return;

        INDEL_PON = Lists.newArrayList();

        // load indel PON
        String refFile = version.is37() ? "/pon/indels_v37.csv" : "/pon/indels_v38.csv";

        final Stream<String> ponLines = new BufferedReader(new InputStreamReader(
                ReferenceData.class.getResourceAsStream(refFile)))
                .lines();

        ponLines.map(Indel::fromString).forEach(INDEL_PON::add);
    }

    public CohortFrequency getAlleleFrequencies() { return mAlleleFrequencies; }

    public HlaAllele findAllele(final String alleleStr, boolean isFourDigit)
    {
        return isFourDigit ? mAlleleCache.findFourDigitAllele(alleleStr) : mAlleleCache.findAllele(alleleStr);
    }

    // HLA gene convenience methods
    public static List<Integer> getAminoAcidExonBoundaries(final HlaGene gene)
    {
        return GENE_CACHE.AminoAcidExonBoundaries.get(gene);
    }

    public static List<Integer> getNucleotideExonBoundaries(final HlaGene gene)
    {
        return GENE_CACHE.NucleotideExonBoundaries.get(gene);
    }

    private void populateAminoAcidSequenceLookup()
    {
        for(HlaSequenceLoci seq : AminoAcidSequences)
        {
            HlaAllele allele = seq.Allele;
            if(!allele.equals(allele.asFourDigit()))
                throw new RuntimeException(format("allele(%s) is not four-digit", allele));

            if(allele.Gene.isPseudo())
                continue;

            ExonSequencesLookup.computeIfAbsent(allele, k -> HlaExonSequences.create(GENE_CACHE.AminoAcidExonBoundaries, seq));
        }
    }

    public boolean load()
    {
        if(mResourceDir.isEmpty()) // a condition for unit testing, otherwise is checked by config loading validation
            return true;

        String nucleotideFilename = mResourceDir + NUC_REF_FILE;

        LL_LOGGER.info("reading nucleotide file: {}", nucleotideFilename);

        if(!loadSequenceFile(nucleotideFilename, NucleotideSequences, false))
            return false;

        String aminoAcidFilename = mResourceDir + AA_REF_FILE;

        LL_LOGGER.info("reading protein file: {}", aminoAcidFilename);

        if(!loadSequenceFile(aminoAcidFilename, AminoAcidSequences, true))
            return false;

        populateAminoAcidSequenceLookup();

        Set<HlaAllele> allelesWithFreqs = Sets.newHashSet(mAlleleFrequencies.getAlleleFrequencies().keySet());
        for(HlaAllele allele : allelesWithFreqs)
        {
            if(!allele.equals(allele.asFourDigit()))
                throw new RuntimeException(format("allele(%s) is not four-digit", allele));

            if(ExonSequencesLookup.containsKey(allele))
                continue;

            LL_LOGGER.warn("allele({}) with cohort frequency has no loaded sequences, dropping from allele frequencies", allele.toString());
            mAlleleFrequencies.getAlleleFrequencies().remove(allele);
        }

        // load and register configured and known alleles
        mAlleleCache.rebuildProteinAlleles(mConfig.ActualAlleles);
        mAlleleCache.rebuildProteinAlleles(mConfig.RestrictedAlleles);

        loadCommonAlleles();

        loadStopLossRecoveryAllele();

        if(HlaYNucleotideSequences != null)
        {
            HlaYNucleotideSequences.addAll(NucleotideSequences.stream().filter(x -> x.Allele.Gene == HLA_Y).toList());
            HlaYNucleotideSequences.forEach(NucleotideSequences::remove);
        }

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
        if(HlaYAminoAcidSequences == null)
            return;

        // construct the AA allele sequences for HLA-Y from the nucleotides if it wasn't loaded
        HlaYAminoAcidSequences.addAll(AminoAcidSequences.stream().filter(x -> x.Allele.Gene == HLA_Y).toList());

        if(!HlaYAminoAcidSequences.isEmpty())
        {
            HlaYAminoAcidSequences.forEach(AminoAcidSequences::remove);
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
            if(sequence.Allele.Gene.isPseudo())
                continue;

            List<Integer> exonBoundaries = getAminoAcidExonBoundaries(sequence.Allele.Gene);
            sequence.setExonBoundaryWildcards(exonBoundaries);
        }

        for(HlaSequenceLoci sequence : NucleotideSequences)
        {
            if(sequence.Allele.Gene.isPseudo())
                continue;

            List<Integer> exonBoundaries = getNucleotideExonBoundaries(sequence.Allele.Gene);
            sequence.setExonBoundaryWildcards(exonBoundaries);
        }
    }

    private void loadCommonAlleles()
    {
        mAlleleFrequencies.getAlleleFrequencies().entrySet().stream()
                .filter(x -> x.getValue() >= COMMON_ALLELES_FREQ_CUTOFF)
                .map(x -> mAlleleCache.requestFourDigit(x.getKey().toString()))
                .forEach(CommonAlleles::add);

        if(!CommonAlleles.isEmpty())
        {
            LL_LOGGER.info("loaded {} common alleles", CommonAlleles.size());
        }
    }

    private void loadStopLossRecoveryAllele()
    {
        // TODO: load from resource file
        if(CURRENT_GENES.coversMhcClass1())
        {
            STOP_LOSS_ON_C_INDEL = mConfig.RefGenVersion.is38() ?
                    new Indel(HLA_CHR, 31269338, "CN", "C") : new Indel(HLA_CHR, 31237115, "CN", "C");

            KnownStopLossIndelAlleles.put(STOP_LOSS_ON_C_INDEL, mAlleleCache.requestFourDigit(STOP_LOSS_ON_C_ALLELE));
        }
    }

    private boolean excludeAllele(final HlaAllele allele)
    {
        if(allele.Gene == HLA_H)
            return true;

        final HlaAllele allele4d = allele.asFourDigit();

        if(EXCLUDED_ALLELES.stream().anyMatch(allele4d::matches))
            return true;

        if(mConfig == null)
            return false;

        if(!mConfig.RestrictedAlleles.isEmpty())
        {
            if(mConfig.ActualAlleles.stream().anyMatch(x -> x.matches(allele4d)))
                return false;

            return mConfig.RestrictedAlleles.stream().noneMatch(x -> x.matches(allele4d));
        }

        return false;
    }

    public static void populateHlaTranscripts(
            final Map<HlaGene, TranscriptData> hlaTranscriptMap, final RefGenomeVersion refGenomeVersion, final GeneSelector genes)
    {
        hlaTranscriptMap.clear();
        hlaTranscriptMap.putAll(loadHlaTranscripts(refGenomeVersion, genes));
    }

    public static Map<HlaGene, TranscriptData> loadHlaTranscripts(
            final RefGenomeVersion refGenomeVersion, @Nullable final GeneSelector genes)
    {
        Map<HlaGene, TranscriptData> hlaTranscriptMap = Maps.newHashMap();

        String transcriptsFile = refGenomeVersion.is37() ? "/transcripts/hla_transcripts_v37.csv" : "/transcripts/hla_transcripts_v38.csv";

        final List<String> hlaTranscriptData = new BufferedReader(new InputStreamReader(
                ReferenceData.class.getResourceAsStream(transcriptsFile))).lines().collect(Collectors.toList());

        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(hlaTranscriptData.get(0), ENSEMBL_DELIM);
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
        TranscriptData currentTrans;
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
                        codingStart, codingEnd, items[biotypeIndex], null);

                HlaGene gene = HlaGene.fromString(geneName);
                if(genes == null || genes.contains(gene))
                    hlaTranscriptMap.put(gene, currentTrans);

                exonDataList = currentTrans.exons();
            }

            ExonData exonData = new ExonData(
                    transId, Integer.parseInt(items[exonStartIndex]), Integer.parseInt(items[exonEndIndex]),
                    Integer.parseInt(items[exonRankIndex]), Integer.parseInt(items[exonPhaseIndex]), Integer.parseInt(items[exonEndPhaseIndex]));

            exonDataList.add(exonData);
        }

        return hlaTranscriptMap;
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
        catch(IOException e)
        {
            LL_LOGGER.error("failed to load ref sequence data from file({}): {}", filename, e.toString());
            return false;
        }
    }

    public void loadSequenceFile(final Iterable<String> fileContents, final Collection<HlaSequenceLoci> sequenceData, boolean isProteinFile)
    {
        for(String line : fileContents)
        {
            String[] items = line.split(",");

            if(items.length != 2)
                return;

            String alleleStr = items[0];

            HlaAllele allele = isProteinFile ? mAlleleCache.requestFourDigit(alleleStr) : mAlleleCache.request(alleleStr);
            if(!CURRENT_GENES.allGenes().contains(allele.Gene))
                continue;

            boolean isDefaultTemplate = isProteinFile && mDeflatedSequenceTemplate == null && allele.matches(DEFLATE_TEMPLATE);
            boolean excludeAllele = excludeAllele(allele);

            if(!isDefaultTemplate && excludeAllele)
                continue;

            String sequenceStr = items[1];
            HlaSequenceLoci newSequence = HlaSequenceFile.createFromReference(allele, sequenceStr, isProteinFile);

            // restrict to first three exons on HLA-DRB1
            if(allele.Gene == HlaGene.HLA_DRB1)
            {
                int nucCount = GENE_CACHE.NucleotideLengths.get(HlaGene.HLA_DRB1);
                int acidCount = nucCount / 3;
                int trimLength = isProteinFile ? acidCount : nucCount;
                List<String> trimmedSequence = Lists.newArrayList(
                        newSequence.getSequences().subList(0, min(trimLength, newSequence.getSequences().size())));
                newSequence = new HlaSequenceLoci(allele, trimmedSequence);
            }

            if(isDefaultTemplate)
                mDeflatedSequenceTemplate = newSequence;

            if(!excludeAllele)
                sequenceData.add(newSequence);
        }
    }
}
