package com.hartwig.hmftools.sage;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.hotspot.ModifiableVariantHotspotEvidence;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageSamConsumer implements Consumer<SAMRecord> {

    /*
    bgzip colo829.sage.vcf
    bcftools index colo829.sage.vcf.gz
    bcftools annotate -a /Users/jon/hmf/resources/GERMLINE_PON.vcf.gz -c GERMLINE_PON_COUNT colo829.sage.vcf.gz -O z -o colo829.sage.germline.vcf.gz
    bcftools annotate -a /Users/jon/hmf/resources/SOMATIC_PON.vcf.gz -c SOMATIC_PON_COUNT colo829.sage.germline.vcf.gz -O z -o colo829.sage.pon.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/out_150_hg19.mappability.bed.gz -h /Users/jon/hmf/resources/mappability.hdr -c CHROM,FROM,TO,-,MAPPABILITY colo829.sage.pon.vcf.gz -O z -o colo829.sage.ann.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -c SOMATIC colo829.sage.ann.vcf.gz -O z -o colo829.sage.final.vcf.gz



    */


    private static final String DISTANCE_FROM_REF_TAG = "NM";
    private final GenomeRegion bounds;
    private final IndexedFastaSequenceFile refGenome;
    private final AtomicInteger count = new AtomicInteger();
    private final ImmutableVariantHotspotImpl.Builder keyBuilder;

    private final Map<VariantHotspot, ModifiableVariantHotspotEvidence> evidenceMap = new HashMap<>();
    private final Map<VariantHotspot, Integer> refSupportMap = new HashMap<>();
    private final Map<VariantHotspot, Integer> refQualityMap = new HashMap<>();

    public SageSamConsumer(@NotNull final GenomeRegion bounds, final IndexedFastaSequenceFile refGenome) {
        this.bounds = bounds;
        this.refGenome = refGenome;
        keyBuilder = ImmutableVariantHotspotImpl.builder().chromosome(bounds.chromosome()).ref("N").alt("N");
    }

    @NotNull
    public List<ModifiableVariantHotspotEvidence> evidence() {

        List<ModifiableVariantHotspotEvidence> result = Lists.newArrayList();
        for (Map.Entry<VariantHotspot, ModifiableVariantHotspotEvidence> entry : evidenceMap.entrySet()) {
            VariantHotspot key = entry.getKey();
            ModifiableVariantHotspotEvidence value = entry.getValue();

            VariantHotspot refKey = ImmutableVariantHotspotImpl.builder().from(key).ref("N").alt("N").build();

            Integer refQuality = refQualityMap.get(refKey);
            if (refQuality != null) {
                value.setIndelSupport(refQuality);
            }

            Integer refSupport = refSupportMap.get(refKey);
            if (refSupport != null) {
                value.setRefSupport(refSupport);
                value.setReadDepth(refSupport + value.altSupport());
            }

            result.add(value);
        }

        Collections.sort(result);
        return result;
    }

    @Override
    public void accept(final SAMRecord record) {

        if (inBounds(record)) {

            List<AlignmentBlock> alignmentBlocks = record.getAlignmentBlocks();
            for (AlignmentBlock alignmentBlock : alignmentBlocks) {

                long refStart = alignmentBlock.getReferenceStart();
                long refEnd = refStart + alignmentBlock.getLength() - 1;
                byte[] refBases = refGenome.getSubsequenceAt(record.getContig(), refStart, refEnd).getBases();

                int readStart = alignmentBlock.getReadStart() - 1;

                String refString = refGenome.getSubsequenceAt(record.getContig(), refStart, refEnd).getBaseString();
                String readString = record.getReadString().substring(readStart, readStart + alignmentBlock.getLength());

                for (int i = 0; i < alignmentBlock.getLength(); i++) {

                    long refPosition = refStart + i;
                    if (!inBounds(refPosition)) {
                        continue;
                    }


                    byte readByte = record.getReadBases()[i + readStart];
                    byte refByte = refBases[i];
                    final int baseQuality = record.getBaseQualities()[i + readStart];

                    if (readByte != refByte) {
                        final VariantHotspot key = keyBuilder.ref(String.valueOf((char) refByte))
                                .alt(String.valueOf((char) readByte))
                                .position(refPosition)
                                .build();

                        ModifiableVariantHotspotEvidence evidence = evidenceMap.computeIfAbsent(key, hotspot -> createEvidence(key));
                        evidence.setAltSupport(evidence.altSupport() + 1);
                        evidence.setAltQuality(evidence.altQuality() + baseQuality);

                        count.incrementAndGet();
                    } else {
//                        final VariantHotspot key = keyBuilder.position(refPosition).ref("N").alt("N").build();
//                        refSupportMap.compute(key, (hotspot, integer) -> integer == null ? 1 : integer + 1);
//                        refQualityMap.compute(key, (hotspot, integer) -> integer == null ? baseQuality : integer + baseQuality);
                    }
                }

            }

        }
    }

    private ModifiableVariantHotspotEvidence createEvidence(VariantHotspot key) {
        return ModifiableVariantHotspotEvidence.create()
                .from(key)
                .setAltSupport(0)
                .setRefSupport(0)
                .setReadDepth(0)
                .setAltQuality(0)
                .setIndelSupport(0);
    }

    public int count() {
        return evidenceMap.size();
    }

    private boolean inBounds(final SAMRecord record) {
        return record.getEnd() >= bounds.start() && record.getStart() <= bounds.end();
    }

    private boolean inBounds(final long position) {
        return position >= bounds.start() && position <= bounds.end();
    }

    public static boolean containsVariant(@NotNull final SAMRecord record) {

        @Nullable
        Object distanceFromReference = record.getAttribute(DISTANCE_FROM_REF_TAG);
        if (distanceFromReference != null) {
            return (int) distanceFromReference != 0;
        }

        return variantInCigar(record.getCigar());
    }

    private static boolean variantInCigar(@NotNull final Cigar cigar) {

        for (final CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator().isIndel() || cigarElement.getOperator().equals(CigarOperator.X)) {
                return true;
            }
        }

        return false;
    }

    private static boolean indelInCigar(@NotNull final Cigar cigar) {

        for (final CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator().isIndel()) {
                return true;
            }
        }

        return false;
    }

    private static boolean clipped(@NotNull final Cigar cigar) {

        for (final CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator().isClipping()) {
                return true;
            }
        }

        return false;
    }

}
