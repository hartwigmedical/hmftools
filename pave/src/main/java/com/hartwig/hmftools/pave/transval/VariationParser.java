package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class VariationParser
{
    static final String HGVS_FORMAT_REQUIRED = "Required format is GENE:p.XnY where X and Y are amino acids and n is an integer.";
    private final EnsemblDataCache ensemblCache;
    private final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap = new HashMap<>();

    public VariationParser(final File ensemblDataDir)
    {
        this.ensemblCache = new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), RefGenomeVersion.V38);
        ensemblCache.setRequiredData(true, true, true, true);
        ensemblCache.load(false);
        ensemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, transcriptAminoAcidsMap, List.of(), false);
    }

    public SingleAminoAcidVariant parse(String input)
    {
        String[] geneVar = input.split(":");
        if (geneVar.length != 2)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        String description = geneVar[1];
        String[] descriptionParts = description.split("\\.");
        if (descriptionParts.length != 2)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        if (!descriptionParts[0].equals("p"))
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }

        Pattern variationPattern = Pattern.compile("^([a-zA-Z]*)(\\d+)([a-zA-Z]+)$");
        final Matcher matcher = variationPattern.matcher(descriptionParts[1]);
        boolean matches = matcher.find();
        if (!matches)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }

        String reference = geneVar[0];
        GeneData geneData = ensemblCache.getGeneDataByName(reference);
        if (geneData == null)
        {
            geneData = ensemblCache.getGeneDataById(reference);
            if (geneData == null)
            {
                throw new IllegalArgumentException(reference + " is not a known gene");
            }
        }
        String referenceAminoAcid = matcher.group(1);
        if (!referenceAminoAcid.isBlank() && isNotValidAminoAcidIdentifier(referenceAminoAcid))
        {
            throw new IllegalArgumentException(referenceAminoAcid + " does not represent an amino acid");
        }
        int position = Integer.parseInt(matcher.group(2));
        String variantAminoAcid = matcher.group(3);
        if (isNotValidAminoAcidIdentifier(variantAminoAcid))
        {
            throw new IllegalArgumentException(variantAminoAcid + " does not represent an amino acid");
        }
        String singleLetterName = AminoAcids.forceSingleLetterProteinAnnotation(variantAminoAcid);

        List<TranscriptData> transcripts = ensemblCache.getTranscripts(geneData.GeneId);
        TranscriptData canonicalTranscript = transcripts.stream()
                .filter(transcriptData -> transcriptData.IsCanonical)
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("No canonical transcript found for " + reference));

        TranscriptAminoAcids aminoAcidsSequence = transcriptAminoAcidsMap.get(canonicalTranscript.TransName);
        if (aminoAcidsSequence == null)
        {
            throw new IllegalStateException("No amino acid sequence found for " + canonicalTranscript.TransName);
        }

        return new SingleAminoAcidVariant(geneData, canonicalTranscript, aminoAcidsSequence, position, singleLetterName);
    }

    private boolean isNotValidAminoAcidIdentifier(String s)
    {
        if (AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s))
        {
            return false;
        }
        if (AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.containsKey(s))
        {
            return false;
        }
        if ("Z".equals(s) || "Glx".equals(s)) // means 'Glutamine or Glutamic Acid'
        {
            return false;
        }
        // B and Asx mean 'Asparagine or Aspartic Acid'
        return !"B".equals(s) && !"Asx".equals(s);
    }
}
