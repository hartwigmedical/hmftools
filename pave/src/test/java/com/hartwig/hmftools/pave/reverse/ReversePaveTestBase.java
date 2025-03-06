package com.hartwig.hmftools.pave.reverse;

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

class ReversePaveTestBase
{
    public final File ensemblDataDir;
    public final ReversePave reversePave;
    public final RefGenomeInterface genome = new TinyGenome();

    public ReversePaveTestBase()
    {
        URL resourceUrl = getClass().getClassLoader().getResource("ensembl_mini");
        try
        {
            assert resourceUrl != null;
            ensemblDataDir = new File(resourceUrl.toURI());
        }
        catch(URISyntaxException e)
        {
            throw new RuntimeException(e);
        }
        reversePave = new ReversePave(ensemblDataDir, RefGenomeVersion.V38, genome);
    }

    @NotNull
    protected static Set<String> css(final String expectedSeparatedByCommas)
    {
        return Arrays.stream(expectedSeparatedByCommas.split(",")).map(String::trim).collect(Collectors.toSet());
    }

    public SingleAminoAcidVariant saav(String definition)
    {
        return (SingleAminoAcidVariant) reversePave.variationParser().parse(definition);
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
        return reversePave.mEnsemblCache.getTranscriptData(geneId, transcriptId);
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

    CodonWithinExons sec(int position, String bases, boolean isForwardStrand)
    {
        return new CodonWithinExons(bs(position, bases, isForwardStrand));
    }

    BaseSequence bs(int position, String bases, boolean isForwardStrand)
    {
        return new BaseSequence(position, bases, isForwardStrand);
    }
}
