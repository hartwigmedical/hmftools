package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class Transval
{
    public final RefGenomeInterface mRefGenomeVersion;
    public final EnsemblDataCache mEnsemblCache;
    private final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap = new HashMap<>();

    public Transval(final File ensemblDataDir, final RefGenomeInterface refGenomeVersion)
    {
        this.mEnsemblCache = new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), RefGenomeVersion.V38);
        this.mRefGenomeVersion = refGenomeVersion;
        mEnsemblCache.setRequiredData(true, true, true, true);
        mEnsemblCache.load(false);
        mEnsemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, transcriptAminoAcidsMap, List.of(), false);
    }

    public VariationParser variationParser()
    {
        return new VariationParser(mEnsemblCache, transcriptAminoAcidsMap);
    }

    public TransvalSNV calculateSNV(String proteinVariant)
    {
        SingleAminoAcidVariant variant = variationParser().parse(proteinVariant);
        String referenceCodon = variant.referenceCodon(mRefGenomeVersion);
        CodonVariantsCalculator change = new CodonVariantsCalculator(variant.referenceAminoAcid(), variant.variantAminoAcid());
        SortedSet<CodonVariant> codonChanges = change.possibleCodonsForVariant(referenceCodon);
        if(codonChanges.isEmpty())
        {
            return null;
        }
        CodonVariant codonVariant = codonChanges.first();
        if(codonVariant.editDistance() != 1)
        {
            return null; // todo deal with this later
        }
        int positionInCodonOfChange = codonVariant.positionOfFirstDifference();
        int positionInChromosomeOfChange = variant.regionsDefiningCodon().translateCodonPosition(positionInCodonOfChange);
        String variantNucleotide = codonVariant.alternateCodon.charAt(positionInCodonOfChange) + "";
        if(!variant.Gene.forwardStrand())
        {
            variantNucleotide = Nucleotides.reverseComplementBases(variantNucleotide);
        }
        String referenceNucleotide = mRefGenomeVersion.getBaseString(variant.Gene.Chromosome, positionInChromosomeOfChange, positionInChromosomeOfChange);
        List<String> alternateNucleotides = new ArrayList<>();
        codonChanges.forEach(cv -> alternateNucleotides.add(cv.alternateCodon));
        return new TransvalSNV(
                variant.Transcript.TransName,
                variant.Gene.Chromosome,
                positionInChromosomeOfChange,
                !variant.codonIsInSingleExon(),
                referenceNucleotide,
                variantNucleotide,
                referenceCodon,
                alternateNucleotides
        );
    }
}
