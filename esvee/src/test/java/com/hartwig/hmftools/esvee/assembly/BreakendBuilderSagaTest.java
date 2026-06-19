package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarAlignedLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignData;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.Breakend;
import com.hartwig.hmftools.esvee.assembly.alignment.BreakendBuilder;
import com.hartwig.hmftools.esvee.assembly.alignment.BreakendSegment;
import com.hartwig.hmftools.esvee.common.saga.SagaAlignment;
import com.hartwig.hmftools.esvee.common.saga.SagaAssembly;
import com.hartwig.hmftools.esvee.common.saga.SagaBreakend;
import com.hartwig.hmftools.esvee.common.saga.SagaVariant;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Test;

import htsjdk.samtools.Cigar;

public class BreakendBuilderSagaTest
{
    private static final String CHROMOSOME = "chr1";
    private static final int ASSEMBLY_LENGTH = 100;
    private static final int SAGA_REF_CONTEXT = 100;
    private static final int SAGA_START = 25;

    private final MockRefGenome mRefGenome;

    private final AssemblyAlignment mAssemblyAlignment;

    public BreakendBuilderSagaTest()
    {
        mRefGenome = new MockRefGenome(true);

        // Since the SAGA variant info is used, nothing in the AssemblyAlignment is used in the code under test. This is mostly dummy data.
        mAssemblyAlignment =
                new AssemblyAlignment(createAssembly(CHROMOSOME, 1000, FORWARD, "A".repeat(ASSEMBLY_LENGTH), ASSEMBLY_LENGTH / 2));
    }

    private void setUpRefGenomeForDel(int lowerBreakend, int upperBreakend, int homologyLength)
    {
        int delLength = upperBreakend - lowerBreakend - 1;
        assertTrue(homologyLength <= delLength);
        String refStart = "A".repeat(lowerBreakend);
        String refDel = "C".repeat(delLength);
        String refHom = "C".repeat(homologyLength);
        String refEnd = "T".repeat(delLength - homologyLength);
        String ref = refStart + refDel + refHom + refEnd;
        assertFalse(mRefGenome.RefGenomeMap.containsKey(CHROMOSOME));
        mRefGenome.RefGenomeMap.put(CHROMOSOME, ref);
    }

    private void setUpRefGenomeForIns(int lowerBreakend, final String insertSeq, int homologyLength)
    {
        assertTrue(homologyLength <= insertSeq.length());
        String refStart = "A".repeat(lowerBreakend);
        // Ref after the insert is homologous to the insert sequence.
        String refHom = insertSeq.substring(0, homologyLength);
        String refEnd = "T".repeat(insertSeq.length() - homologyLength);
        String ref = refStart + refHom + refEnd;
        assertFalse(mRefGenome.RefGenomeMap.containsKey(CHROMOSOME));
        mRefGenome.RefGenomeMap.put(CHROMOSOME, ref);
    }

    private static SagaAssembly deletionSagaAssembly(int lowerPos, int upperPos)
    {
        SagaVariant variant = new SagaVariant(
                "test-del",
                new SagaBreakend(new BasePosition(CHROMOSOME, lowerPos), FORWARD),
                new SagaBreakend(new BasePosition(CHROMOSOME, upperPos), REVERSE),
                ""
        );
        String refStart = "A".repeat(SAGA_REF_CONTEXT);
        String refEnd = "T".repeat(SAGA_REF_CONTEXT);
        String sequence = refStart + refEnd;
        List<Integer> junctionOffsets = List.of(refStart.length());
        return new SagaAssembly("", variant, junctionOffsets, sequence);
    }

    private static SagaAssembly insertionSagaAssembly(int lowerPos, final String insertSequence)
    {
        int upperPos = lowerPos + 1;
        SagaVariant variant = new SagaVariant(
                "test-ins",
                new SagaBreakend(new BasePosition(CHROMOSOME, lowerPos), FORWARD),
                new SagaBreakend(new BasePosition(CHROMOSOME, upperPos), REVERSE),
                insertSequence
        );
        String refStart = "A".repeat(SAGA_REF_CONTEXT);
        String refEnd = "T".repeat(SAGA_REF_CONTEXT);
        String sequence = refStart + insertSequence + refEnd;
        List<Integer> junctionOffsets = List.of(refStart.length(), sequence.length() - refEnd.length());
        return new SagaAssembly("", variant, junctionOffsets, sequence);
    }

    private SagaAlignment sagaAlignment(final SagaAssembly sagaAssembly, int sagaStart, final String cigarStr)
    {
        Cigar cigar = cigarFromStr(cigarStr);
        int queryLength = mAssemblyAlignment.fullSequenceLength();
        int queryStart = leftClipLength(cigar);
        int queryEnd = queryLength - rightClipLength(cigar);
        int sagaEnd = sagaStart + cigarAlignedLength(cigar);
        BwaMemAlignment raw = new BwaMemAlignment(
                0, 0,   // Doesn't matter
                sagaStart, sagaEnd, queryStart, queryEnd,
                0, 0, 0, 0,    // Doesn't matter
                cigarStr,
                "", "", 0, 0, 0 // Doesn't matter
        );
        SagaAlignment sagaAlignment = new SagaAlignment(raw, cigar, queryLength, sagaAssembly);
        sagaAlignment.validate();
        return sagaAlignment;
    }

    private SagaAlignment sagaAlignment(final SagaAssembly sagaAssembly)
    {
        return sagaAlignment(sagaAssembly, SAGA_START, format("%dM", mAssemblyAlignment.fullSequenceLength()));
    }

    private List<Breakend> formBreakends(final SagaAlignment sagaAlignment)
    {
        AlignData alignData = AlignData.fromSaga(sagaAlignment);
        BreakendBuilder builder = new BreakendBuilder(mRefGenome, mAssemblyAlignment);
        assertTrue(builder.formBreakendsFromSaga(sagaAlignment, alignData));
        return mAssemblyAlignment.breakends();
    }

    @Test
    public void testDeletionNoHomology()
    {
        // 200b DEL with no homology to ref.

        int lowerPos = 1001;
        int upperPos = 1200;

        // TODO: don't need to specify chromosome here?
        setUpRefGenomeForDel(lowerPos, upperPos, 0);

        SagaAssembly sagaAssembly = deletionSagaAssembly(lowerPos, upperPos);
        SagaAlignment sagaAlignment = sagaAlignment(sagaAssembly);

        List<Breakend> breakends = formBreakends(sagaAlignment);
        assertEquals(2, breakends.size());

        Breakend lower = breakends.get(0);
        assertEquals(CHROMOSOME, lower.Chromosome);
        assertEquals(lowerPos, lower.Position);
        assertEquals(FORWARD, lower.Orient);
        assertEquals(DEL, lower.svType());
        assertFalse(lower.isSagaInferred());

        Breakend upper = breakends.get(1);
        assertEquals(CHROMOSOME, upper.Chromosome);
        assertEquals(upperPos, upper.Position);
        assertEquals(REVERSE, upper.Orient);
        assertEquals(DEL, upper.svType());
        assertFalse(upper.isSagaInferred());
    }

    @Test
    public void testDeletionWithHomology()
    {
        // 200b DEL with 10b homology to ref.

        int lowerPos = 1001;
        int upperPos = 1200;
        int homology = 10;

        setUpRefGenomeForDel(lowerPos, upperPos, homology);

        SagaAssembly sagaAssembly = deletionSagaAssembly(lowerPos, upperPos);
        SagaAlignment sagaAlignment = sagaAlignment(sagaAssembly);

        List<Breakend> breakends = formBreakends(sagaAlignment);
        assertEquals(2, breakends.size());

        Breakend lower = breakends.get(0);
        assertEquals(CHROMOSOME, lower.Chromosome);
        assertEquals(lowerPos + homology / 2, lower.Position);
        assertEquals(FORWARD, lower.Orient);
        assertEquals(DEL, lower.svType());
        assertFalse(lower.isSagaInferred());

        Breakend upper = breakends.get(1);
        assertEquals(CHROMOSOME, upper.Chromosome);
        assertEquals(upperPos + homology / 2, upper.Position);
        assertEquals(REVERSE, upper.Orient);
        assertEquals(DEL, upper.svType());
        assertFalse(upper.isSagaInferred());
    }

    @Test
    public void testInsertionNoHomology()
    {
        // 20b INS with no homology to ref.

        int lowerPos = 1001;
        String insertSequence = "ACCGTACCGTGAGAGAGAGA";
        int upperPos = lowerPos + 1;
        int homology = 0;

        setUpRefGenomeForIns(lowerPos, insertSequence, homology);

        SagaAssembly sagaAssembly = insertionSagaAssembly(lowerPos, insertSequence);
        SagaAlignment sagaAlignment = sagaAlignment(sagaAssembly);

        List<Breakend> breakends = formBreakends(sagaAlignment);
        assertEquals(2, breakends.size());

        Breakend lower = breakends.get(0);
        assertEquals(CHROMOSOME, lower.Chromosome);
        assertEquals(lowerPos, lower.Position);
        assertEquals(FORWARD, lower.Orient);
        assertEquals(insertSequence, lower.InsertedBases);
        assertEquals(INS, lower.svType());

        Breakend upper = breakends.get(1);
        assertEquals(CHROMOSOME, upper.Chromosome);
        assertEquals(upperPos, upper.Position);
        assertEquals(REVERSE, upper.Orient);
        assertEquals(insertSequence, lower.InsertedBases);
        assertEquals(INS, lower.svType());
    }

    @Test
    public void testInsertionWithHomology()
    {
        // 20b INS with 10b homology to ref. Gets converted to a DUP.

        int lowerPos = 1001;
        String insertSequence = "ACCGTTATATGAGAGAGAGA";
        int upperPos = lowerPos + 1;
        int homology = 10;

        setUpRefGenomeForIns(lowerPos, insertSequence, homology);

        SagaAssembly sagaAssembly = insertionSagaAssembly(lowerPos, insertSequence);
        SagaAlignment sagaAlignment = sagaAlignment(sagaAssembly);

        List<Breakend> breakends = formBreakends(sagaAlignment);
        assertEquals(2, breakends.size());

        Breakend lower = breakends.get(0);
        assertEquals(CHROMOSOME, lower.Chromosome);
        assertEquals(lowerPos + 1, lower.Position);
        assertEquals(REVERSE, lower.Orient);
        String expectedInsSeq = insertSequence.substring(homology);
        assertEquals(expectedInsSeq, lower.InsertedBases);
        assertEquals(DUP, lower.svType());
        ;

        Breakend upper = breakends.get(1);
        assertEquals(CHROMOSOME, upper.Chromosome);
        assertEquals(upperPos + homology - 1, upper.Position);
        assertEquals(FORWARD, upper.Orient);
        assertEquals(expectedInsSeq, lower.InsertedBases);
        assertEquals(DUP, lower.svType());
    }

    @Test
    public void testSequenceIndicesWithCigarIndel()
    {
        // SAGA INS: chr2:200→201, insertSequence="ACGTACGT" (8 bases).
        // junctionOffsets in the SAGA sequence are [50, 58].
        // The assembly has a 2-base insertion relative to SAGA at position 50 (cigar "50M2I58M"),
        // so the query junction offsets are shifted by +2 vs. a simple 1:1 mapping:
        //   sagaIndex=50 → queryIndex=52
        //   sagaIndex=58 → queryIndex=60
        // The BreakendSegment.indelSeqenceIndices must reflect these query-space positions.

        String insertSequence = "ACGTACGT";
        SagaAssembly sagaAssembly = insertionSagaAssembly(200, insertSequence);
        // "50M2I58M": 50+58=108 ref bases (matches SAGA length 108), query length = 50+2+58=110
        SagaAlignment sagaAlignment = sagaAlignment(sagaAssembly, 0, "50M2I58M");

        List<Breakend> breakends = formBreakends(sagaAlignment);
        assertEquals(2, breakends.size());

        // Both breakends share the same segment; retrieve it from either.
        BreakendSegment segment = breakends.get(0).segments().get(0);
        // SequenceIndex = alignData.sequenceStart() = 0.
        assertEquals(0, segment.SequenceIndex);
        // Query junction indices: 52 and 60 (each offset by +2 from the cigar insertion).
        int[] indelIndices = segment.indelSeqenceIndices();
        assertEquals(52, indelIndices[0]);
        assertEquals(60, indelIndices[1]);
    }

    @Test
    public void testInferredBreakends()
    {
        // SAGA INS: chr2:200→201, insertSequence="ACGTACGT", junctionOffsets=[50, 58].
        // Assembly only aligns to SAGA positions [55, 57): both junctions (50 and 58) are outside
        // this range, so both breakends are inferred from SAGA rather than from the assembly itself.
        //   startInferredLength = sagaStart(55) - junctionOffset(50) = 5 > 0 → lower inferred
        //   endInferredLength   = junctionOffset(58) - sagaEnd(57)   = 1 > 0 → upper inferred

        String insertSequence = "ACGTACGT";
        SagaAssembly sagaAssembly = insertionSagaAssembly(200, insertSequence);
        SagaAlignment sagaAlignment = sagaAlignment(sagaAssembly, 55, "2M");

        List<Breakend> breakends = formBreakends(sagaAlignment);
        assertEquals(2, breakends.size());

        assertTrue(breakends.get(0).isSagaInferred());
        assertTrue(breakends.get(1).isSagaInferred());
    }
}
