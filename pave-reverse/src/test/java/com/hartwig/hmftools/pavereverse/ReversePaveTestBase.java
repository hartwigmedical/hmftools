package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.jetbrains.annotations.NotNull;

import com.hartwig.hmftools.common.utils.EnsemblMini;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcid;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSpecification;
import com.hartwig.hmftools.pavereverse.base.BaseSequence;
import com.hartwig.hmftools.pavereverse.base.CodonWithinExons;
import com.hartwig.hmftools.pavereverse.base.SplitCodonSequence;
import com.hartwig.hmftools.pavereverse.protein.SingleAminoAcidVariant;

public class ReversePaveTestBase
{
    public final File ensemblDataDir;
    public final ReversePave reversePave;
    public final RefGenomeInterface genome = new TinyGenome();
    protected final String braf = "BRAF";
    protected final String brafCanonical = "ENST00000646891";
    protected final String zyx = "ZYX";
    protected final String zyxCanonical = "ENST00000322764";
    protected final String tatdn2 = "TATDN2";
    protected final String tatdn2Canonical = "ENST00000448281";

    public ReversePaveTestBase()
    {
        ensemblDataDir = EnsemblMini.ensemblMiniDataDir();
        reversePave = new ReversePave(ensemblDataDir, RefGenomeVersion.V38, genome);
    }

    @NotNull
    protected static Set<String> css(final String expectedSeparatedByCommas)
    {
        return Arrays.stream(expectedSeparatedByCommas.split(",")).map(String::trim).collect(Collectors.toSet());
    }

    public SingleAminoAcidVariant saav(String definition)
    {
        return (SingleAminoAcidVariant) reversePave.proteinVariantParser().parse(definition);
    }

    protected SplitCodonSequence seq(String left, String right)
    {
        return new SplitCodonSequence(left, right, 100);
    }

    protected BaseSequenceChange basesChange(String ref, String alt, String chr, int position)
    {
        return new BaseSequenceChange(ref, alt, chr, position);
    }

    protected AminoAcid aa(String s)
    {
        return new AminoAcid(s);
    }

    protected TranscriptData transcript(String geneId, String transcriptId)
    {
        return reversePave.EnsemblCache.getTranscriptData(geneId, transcriptId);
    }

    protected AminoAcidSpecification aas(int position, String symbol)
    {
        final AminoAcid aa = symbol == null ? null : aa(symbol);
        return new AminoAcidSpecification(position, aa);
    }

    protected AminoAcidSequence aaSeq(String sequence)
    {
        return AminoAcidSequence.parse(sequence);
    }

    protected void check(BaseSequenceChange change, String ref, String alt, String chr, int position)
    {
        assertEquals(ref, change.Ref);
        assertEquals(alt, change.Alt);
        assertEquals(chr, change.Chromosome);
        assertEquals(position, change.Position);
    }

    protected void checkSingleChange(BaseSequenceVariants variant, String ref, String alt, String chr, int position)
    {
        Set<BaseSequenceChange> hotspots = variant.changes();
        assertEquals(1, hotspots.size());
        assertTrue(hotspots.contains(basesChange(ref, alt, chr, position)));
    }

    protected void checkChanges(BaseSequenceVariants variant, BaseSequenceChange... expected)
    {
        Set<BaseSequenceChange> hotspots = variant.changes();
        assertEquals(expected.length, hotspots.size());
        for(BaseSequenceChange hotspot : expected)
        {
            assertTrue(hotspots.contains(hotspot));
        }
    }

    protected CodonWithinExons sec(int position, String bases, boolean isForwardStrand)
    {
        return new CodonWithinExons(bs(position, bases, isForwardStrand));
    }

    public BaseSequence bs(int position, String bases, boolean isForwardStrand)
    {
        return new BaseSequence(position, bases, isForwardStrand);
    }
}
