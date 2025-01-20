package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.lang3.tuple.Pair;

public class Transval
{
    public final RefGenomeInterface mRefGenomeVersion;
    public final EnsemblDataCache mEnsemblCache;
    private final Map<String, TranscriptAminoAcids> mTranscriptAminoAcidsMap = new HashMap<>();

    public Transval(final File ensemblDataDir, final RefGenomeInterface refGenomeVersion)
    {
        this.mEnsemblCache = new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), RefGenomeVersion.V38);
        this.mRefGenomeVersion = refGenomeVersion;
        mEnsemblCache.setRequiredData(true, true, true, true);
        mEnsemblCache.load(false);
        mEnsemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTranscriptAminoAcidsMap, List.of(), false);
    }

    public VariationParser variationParser()
    {
        return new VariationParser(mEnsemblCache, mTranscriptAminoAcidsMap);
    }

    public TransvalSnvMnv calculateVariant(String proteinVariant)
    {
        SingleAminoAcidVariant variant = variationParser().parse(proteinVariant);
        String referenceCodon = variant.referenceCodon(mRefGenomeVersion);
        CodonVariantsCalculator change = new CodonVariantsCalculator(variant.referenceAminoAcids(), variant.altValue());
        SortedSet<CodonVariant> codonChanges = change.possibleCodonsForVariant(referenceCodon);
        if(codonChanges.isEmpty())
        {
            return null;
        }
        CodonVariant codonVariant = codonChanges.first();
        int positionInCodonOfChange = codonVariant.positionOfFirstDifference();
        int positionInChromosomeOfChange = variant.regionsDefiningCodon().translateCodonPosition(positionInCodonOfChange);
        Pair<String,String> nucleotideDifferences = codonVariant.differenceStrings();
        if(!variant.Gene.forwardStrand())
        {
            nucleotideDifferences = Pair.of(reverseComplementBases(nucleotideDifferences.getLeft()), reverseComplementBases(nucleotideDifferences.getRight()));
        }
        List<String> alternateNucleotides = new ArrayList<>();
        codonChanges.forEach(cv -> alternateNucleotides.add(cv.AlternateCodon));
        return new TransvalSnvMnv(
                variant.Transcript.TransName,
                variant.Gene.Chromosome,
                positionInChromosomeOfChange,
                !variant.codonIsInSingleExon(),
                nucleotideDifferences.getLeft(),
                nucleotideDifferences.getRight(),
                referenceCodon,
                alternateNucleotides
        );
    }
}
