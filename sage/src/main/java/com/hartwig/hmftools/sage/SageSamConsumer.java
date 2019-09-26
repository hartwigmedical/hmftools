package com.hartwig.hmftools.sage;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hotspot.ModifiableVariantHotspotEvidence;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.sage.count.BaseDetails;
import com.hartwig.hmftools.sage.count.EvictingLinkedMap;
import com.hartwig.hmftools.sage.count.ReadContext;

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

    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA -c FILTER colo829.sage.map.vcf.gz -O z -o colo829.sage.pre.vcf.gz
    bcftools index colo829.sage.pre.vcf.gz

    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA -c FILTER colo829.sage.pre.vcf.gz -O z -o colo829.sage.final.vcf.gz
    */

    private static final String DISTANCE_FROM_REF_TAG = "NM";
    private final int minQuality;

    private final GenomeRegion bounds;
    private final IndexedFastaSequenceFile refGenome;

    private final EvictingLinkedMap<Long, BaseDetails> baseMap;
    private final List<BaseDetails> baseList = Lists.newArrayList();

    public SageSamConsumer(final int minQuality, @NotNull final GenomeRegion bounds, final IndexedFastaSequenceFile refGenome) {
        this(minQuality, bounds, refGenome, Sets.newHashSet());

    }

    public SageSamConsumer(final int minQuality, @NotNull final GenomeRegion bounds, final IndexedFastaSequenceFile refGenome,
            Set<Long> hotspots) {
        this.bounds = bounds;
        this.refGenome = refGenome;
        this.minQuality = minQuality;

        final BiConsumer<Long, BaseDetails> baseDetailsHandler = (position, baseDetails) -> {
            if (!baseDetails.isEmpty() || hotspots.contains(baseDetails.position())) {
                baseList.add(baseDetails);
            }
        };
        baseMap = new EvictingLinkedMap<>(baseDetailsHandler);
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

            if (record.getMappingQuality() >= minQuality) {
                record.getAlignmentBlocks().forEach(x -> processPrimeAlignment(record, x));
            } else {
                record.getAlignmentBlocks().forEach(x -> processSubprimeAlignment(record, x));
            }

        }
    }

    private void processSubprimeAlignment(@NotNull final SAMRecord record, @NotNull final AlignmentBlock alignmentBlock) {
        long refStart = alignmentBlock.getReferenceStart();
        for (int refBytePosition = 0; refBytePosition < alignmentBlock.getLength(); refBytePosition++) {
            long position = refStart + refBytePosition;
            baseDetails(record.getContig(), position).incrementSubprimeReadDepth();
        }
    }

    private void processPrimeAlignment(@NotNull final SAMRecord record, @NotNull final AlignmentBlock alignmentBlock) {

        long refStart = alignmentBlock.getReferenceStart();
        long refEnd = refStart + alignmentBlock.getLength() - 1;
        byte[] refBytes = refGenome.getSubsequenceAt(record.getContig(), refStart, refEnd).getBases();

        int readStart = alignmentBlock.getReadStart() - 1;

        for (int refBytePosition = 0; refBytePosition < alignmentBlock.getLength(); refBytePosition++) {

            long position = refStart + refBytePosition;
            int readBytePosition = readStart + refBytePosition;

            if (!inBounds(position)) {
                continue;
            }

            final byte refByte = refBytes[refBytePosition];
            final byte readByte = record.getReadBases()[readBytePosition];
            final int baseQuality = record.getBaseQualities()[readBytePosition];

            final BaseDetails baseDetails = baseDetails(record.getContig(), position);

            baseDetails.incrementReadDepth();

            if (readByte != refByte) {

                final String ref = String.valueOf((char) refByte);
                final String alt = String.valueOf((char) readByte);

                long distanceFromAlignmentStart = position - record.getAlignmentStart();
                long distanceFromAlignmentEnd = record.getAlignmentEnd() - position;
                int minDistanceFromAlignment = (int) Math.min(distanceFromAlignmentStart, distanceFromAlignmentEnd);

                final ModifiableVariantHotspotEvidence evidence = baseDetails.selectOrCreate(ref, alt);
                evidence.setAltSupport(evidence.altSupport() + 1);
                evidence.setAltQuality(evidence.altQuality() + baseQuality);
                evidence.setAltMapQuality(evidence.altMapQuality() + record.getMappingQuality());
                evidence.setAltMinQuality(evidence.altMinQuality() + Math.min(record.getMappingQuality(), baseQuality));
                evidence.setAltDistanceFromRecordStart(evidence.altDistanceFromRecordStart() + readBytePosition);
                evidence.setAltMinDistanceFromAlignment(evidence.altMinDistanceFromAlignment() + minDistanceFromAlignment);

                baseDetails.addReadContext(new ReadContext(readBytePosition, record.getReadBases()));

            } else {
                baseDetails.incrementRefSupport();
                baseDetails.incrementRefQuality(baseQuality);
            }
        }

    }

    @NotNull
    private BaseDetails baseDetails(@NotNull final String contig, final long position) {
        return baseMap.compute(position, (key, old) -> old == null ? new BaseDetails(contig, position) : old);
    }

    private int mismatchedBases(int refStart, int readStart, byte[] refBases, byte[] readBases) {
        assert (refBases[refStart] != readBases[readStart]);


        //TODO: This feels a little dodgy
        // TODO: Fix up quality for MNVs && Increment read depth of subsequent bases
        //                        int mismatchedBases = mismatchedBases(refBytePosition, readBytePosition, refBytes, record.getReadBases());
        //                        final String ref = new String(Arrays.copyOfRange(refBytes, refBytePosition, refBytePosition + mismatchedBases));
        //                        final String alt =
        //                                new String(Arrays.copyOfRange(record.getReadBases(), readBytePosition, readBytePosition + mismatchedBases));
        //                        refBytePosition += (mismatchedBases - 1);
        //                        for (int i = 1; i < mismatchedBases; i++) {
        //                            baseDetails(record.getContig(), position + i).incrementReadDepth();
        //                        }

        for (int i = 1; ; i++) {
            int refPosition = refStart + i;
            int readPosition = readStart + i;

            if (refPosition >= refBases.length || readPosition >= readBases.length || refBases[refPosition] == readBases[readPosition]) {
                return i;
            }
        }
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
