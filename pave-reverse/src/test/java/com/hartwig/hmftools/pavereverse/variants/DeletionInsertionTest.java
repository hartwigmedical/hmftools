package com.hartwig.hmftools.pavereverse.variants;

import static org.junit.Assert.assertEquals;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;
import com.hartwig.hmftools.pavereverse.base.PaddedExon;

import org.junit.Assert;
import org.junit.Test;

public class DeletionInsertionTest extends VariantTest
{
    @Test
    public void variantSequenceTest()
    {
        AminoAcidSequence replacement = AminoAcidSequence.parse("TKWRF");
        DeletionInsertion di = new DeletionInsertion(gene, transcript, taa, aar, replacement);
        // MAAQVAPAAS -> MA + TKWRF + VAPAAS
        AminoAcidSequence expected = AminoAcidSequence.parse("MATKWRFVAPAAS");
        assertEquals(expected, di.variantSequence());
    }

    @Test
    public void seekResultsInCompanionContextTest()
    {
        AminoAcidSequence replacement = AminoAcidSequence.parse("TKWRF");
        DeletionInsertion di = new DeletionInsertion(gene, transcript, taa, aar, replacement);
        Assert.assertFalse(di.seekResultsInCompanionContext(true));
    }

    @Test
    public void selectChangesToReportTest()
    {
        // Exon0 bases: ATG GCC GCG CAG GTC
        AminoAcidSequence replacement = AminoAcidSequence.parse("AE");
        DeletionInsertion di = new DeletionInsertion(gene, transcript, taa, aar, replacement);
        Set<ChangeResult> changes = new HashSet<>();
        changes.add(new ChangeResult(replacement, exon0Bases, 6, "GCGCAG", "GCTGAA"));
        changes.add(new ChangeResult(replacement, exon0Bases, 6, "GCGCAG", "GCCGAG"));
        changes.add(new ChangeResult(replacement, exon0Bases, 6, "GCGCAG", "GCAGAA"));
        Collection<ChangeResult> reportingChanges = di.selectChangesToReport(changes);
        assertEquals(1, reportingChanges.size());
        assertEquals("GCCGAG", reportingChanges.iterator().next().AltBases);
    }

    @Test
    public void veryLongInsertionTest()
    {
        final String longSeq = "LRSALRSLRSLRSTAALRSTLL";
        AminoAcidSequence replacement = AminoAcidSequence.parse(longSeq);

        PaddedExon exon = new PaddedExon(8,"", "", exon0Bases, 9, "GGATC", "TACG");
        ChangeContext context = new ChangeContext(exon, 9, 10, true, 1);
        DeletionInsertion di = new DeletionInsertion(gene, transcript, taa, aar, replacement);
        Set<ChangeResult> results = di.applyChange(context);
        Assert.assertEquals(1, results.size());
        ChangeResult result = results.iterator().next();
        Assert.assertEquals("MA" + longSeq + "V", result.Acids.sequence());
    }
}
