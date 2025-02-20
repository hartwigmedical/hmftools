package com.hartwig.hmftools.pave.transval;

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

import org.jetbrains.annotations.NotNull;

class TransvalTestBase
{
    public final File ensemblDataDir;
    public final Transval transval;
    public final RefGenomeInterface genome = new TinyGenome();

    public TransvalTestBase()
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
        transval = new Transval(ensemblDataDir, genome);
    }

    @NotNull
    protected static Set<String> css(final String expectedSeparatedByCommas)
    {
        return Arrays.stream(expectedSeparatedByCommas.split(",")).map(String::trim).collect(Collectors.toSet());
    }

    public SingleAminoAcidVariant saav(String definition)
    {
        return transval.variationParser().parseSingleAminoAcidVariant(definition);
    }

    protected SplitCodonSequence seq(String left, String right)
    {
        return new SplitCodonSequence(left, right, 100);
    }

    protected TransvalHotspot hotspot(String ref, String alt, String chr, int position)
    {
        return new TransvalHotspot(ref, alt, chr, position);
    }

    protected AminoAcid aa(String s)
    {
        return new AminoAcid(s);
    }

    protected TranscriptData transcript(String geneId, String transcriptId)
    {
        return transval.mEnsemblCache.getTranscriptData(geneId, transcriptId);
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

    protected void checkSingleHotspot(TransvalVariant variant, String ref, String alt, String chr, int position)
    {
        Set<TransvalHotspot> hotspots = variant.hotspots();
        assertEquals(1, hotspots.size());
        assertTrue(hotspots.contains(hotspot(ref, alt, chr, position)));
    }

    protected void checkHotspots(TransvalVariant variant, TransvalHotspot... expected)
    {
        Set<TransvalHotspot> hotspots = variant.hotspots();
        assertEquals(expected.length, hotspots.size());
        for(TransvalHotspot hotspot : expected)
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
