package com.hartwig.hmftools.sage;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hotspot.ModifiableVariantHotspotEvidence;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.sage.count.BaseDetails;
import com.hartwig.hmftools.sage.count.EvictingLinkedMap;

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
    bcftools index colo829.sage.germline.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/SOMATIC_PON.vcf.gz -c SOMATIC_PON_COUNT colo829.sage.germline.vcf.gz -O z -o colo829.sage.pon.vcf.gz
    bcftools index colo829.sage.pon.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/out_150_hg19.mappability.bed.gz -h /Users/jon/hmf/resources/mappability.hdr -c CHROM,FROM,TO,-,MAPPABILITY colo829.sage.pon.vcf.gz -O z -o colo829.sage.map.vcf.gz
    bcftools index colo829.sage.map.vcf.gz

    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m STRELKA colo829.sage.map.vcf.gz -O z -o colo829.sage.final.vcf.gz


    */

    private static final String DISTANCE_FROM_REF_TAG = "NM";
    private final GenomeRegion bounds;
    private final IndexedFastaSequenceFile refGenome;

    private final Set<Long> hotspots;
    private final EvictingLinkedMap<Long, BaseDetails> baseMap;
    private final List<BaseDetails> baseList = Lists.newArrayList();

    public SageSamConsumer(@NotNull final GenomeRegion bounds, final IndexedFastaSequenceFile refGenome) {
        this(bounds, refGenome, Sets.newHashSet());

    }

    public SageSamConsumer(@NotNull final GenomeRegion bounds, final IndexedFastaSequenceFile refGenome, Set<Long> hotspots) {
        this.bounds = bounds;
        this.refGenome = refGenome;

        BiConsumer<Long, BaseDetails> baseDetailsHandler = (position, baseDetails) -> {
            if (!baseDetails.isEmpty() || hotspots.contains(baseDetails.position())) {
                baseList.add(baseDetails);
            }

        };
        baseMap = new EvictingLinkedMap<>(baseDetailsHandler);
        this.hotspots = hotspots;
    }

    @NotNull
    public List<BaseDetails> bases() {
        baseMap.evictAll();
        Collections.sort(baseList);
        return baseList;
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

                for (int i = 0; i < alignmentBlock.getLength(); i++) {
                    long refPosition = refStart + i;
                    int readPosition = readStart + i;

                    if (!inBounds(refPosition)) {
                        continue;
                    }

                    byte refByte = refBases[i];
                    byte readByte = record.getReadBases()[readPosition];

                    final String ref = String.valueOf((char) refByte);
                    final String alt = String.valueOf((char) readByte);

                    final BaseDetails baseDetails = baseMap.compute(refPosition,
                            (position, old) -> old == null ? new BaseDetails(record.getContig(), position) : old);

                    baseDetails.incrementReadDepth();
                    baseDetails.incrementDistanceFromRecordStart(readPosition);
                    baseDetails.incrementRecordDistances((int)refPosition, record.getAlignmentStart(), record.getAlignmentEnd());

                    final int baseQuality = record.getBaseQualities()[i + readStart];

                    if (readByte != refByte) {
                        final ModifiableVariantHotspotEvidence evidence = baseDetails.selectOrCreate(ref, alt);
                        evidence.setAltSupport(evidence.altSupport() + 1);
                        evidence.setAltQuality(evidence.altQuality() + baseQuality);

                    } else {
                        baseDetails.incrementRefSupport();
                        baseDetails.incrementRefQuality(baseQuality);
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
