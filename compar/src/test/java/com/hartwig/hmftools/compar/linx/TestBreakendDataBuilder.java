package com.hartwig.hmftools.compar.linx;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestBreakendDataBuilder
{
    public String transcriptId = "ENST00000332149";
    public boolean canonical = true;
    public boolean reportedDisruption = true;
    public TranscriptRegionType regionType = TranscriptRegionType.INTRONIC;
    public TranscriptCodingType codingType = TranscriptCodingType.CODING;
    public int nextSpliceExonRank = 2;
    public StructuralVariantType svType = StructuralVariantType.BND;
    public String chromosome = "chr21";
    public int position = 41500000;
    public byte orientation = 0;
    public int[] homologyOffset = {0, 0};
    public String comparisonChromosome = "chr21";
    public int comparisonPosition = 41500000;

    private static final Consumer<TestBreakendDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.transcriptId = "ENST00000297405";
        b.canonical = false;
        b.reportedDisruption = false;
        b.regionType = TranscriptRegionType.EXONIC;
        b.codingType = TranscriptCodingType.NON_CODING;
        b.nextSpliceExonRank = 7;
        b.svType = StructuralVariantType.DUP;
        b.chromosome = "chr1";
        b.position = 10000;
        b.orientation = 1;
        b.homologyOffset = new int[] { -9, 10 };
        b.comparisonChromosome = "chr1";
        b.comparisonPosition = 10000;
    };

    public static final TestComparableItemBuilder<TestBreakendDataBuilder, BreakendData> BUILDER =
            new TestComparableItemBuilder<>(TestBreakendDataBuilder::new, TestBreakendDataBuilder::build, ALTERNATE_INITIALIZER);

    private BreakendData build()
    {
        final LinxBreakend breakend = ImmutableLinxBreakend.builder()
                .transcriptId(transcriptId)
                .canonical(canonical)
                .reportedDisruption(reportedDisruption)
                .regionType(regionType)
                .codingType(codingType)
                .nextSpliceExonRank(nextSpliceExonRank)
                .id(-1)
                .svId(-1)
                .isStart(true)
                .gene("")
                .geneOrientation("")
                .disruptive(true)
                .undisruptedCopyNumber(-1)
                .biotype("")
                .exonUp(-1)
                .exonDown(-1)
                .exonicBasePhase(-1)
                .nextSpliceExonPhase(-1)
                .nextSpliceDistance(-1)
                .totalExonCount(-1)
                .build();
        return new BreakendData(
                breakend, "", svType, chromosome, position, orientation, homologyOffset, comparisonChromosome, comparisonPosition);
    }
}
