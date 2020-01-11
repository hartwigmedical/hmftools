package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.read.ReadContextFactory.createDelContext;
import static com.hartwig.hmftools.sage.read.ReadContextFactory.createInsertContext;
import static com.hartwig.hmftools.sage.read.ReadContextFactory.createMNVContext;
import static com.hartwig.hmftools.sage.read.ReadContextFactory.createSNVContext;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.sam.CigarHandler;
import com.hartwig.hmftools.sage.sam.CigarTraversal;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class RefContextConsumer implements Consumer<SAMRecord> {

    private static final Logger LOGGER = LogManager.getLogger(RefContextConsumer.class);

    /*

    bcftools filter -i 'FORMAT/QUAL[1:0]>99 && FORMAT/AD[0:1] < 4' GIABvsSELFv004.sage.vcf.gz -O z -o GIABvsSELFv004.sage.filtered.vcf.gz
    bcftools filter -i 'FORMAT/QUAL[1:0]>99 && FORMAT/AD[0:1] < 4' COLO829v003.sage.vcf.gz -O z -o COLO829v003.sage.filtered.vcf.gz



    bgzip colo829.sage.vcf
    bcftools index colo829.sage.vcf.gz

    input=$1

    bcftools annotate -a /Users/jon/hmf/resources/GERMLINE_PON.vcf.gz -c GERMLINE_PON_COUNT ${input}.vcf.gz -O z -o ${input}.germline.vcf.gz
    bcftools index ${input}.germline.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/SOMATIC_PON.vcf.gz -c SOMATIC_PON_COUNT ${input}.germline.vcf.gz -O z -o ${input}.pon.vcf.gz
    bcftools index ${input}.pon.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/out_150_hg19.mappability.bed.gz -h /Users/jon/hmf/resources/mappability.hdr -c CHROM,FROM,TO,-,MAPPABILITY ${input}.pon.vcf.gz -O z -o ${input}.map.vcf.gz
    bcftools index ${input}.map.vcf.gz

    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA -c FILTER colo829.sage.map.vcf.gz -O z -o colo829.sage.pre.vcf.gz
    bcftools index colo829.sage.pre.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA -c FILTER colo829.sage.pre.vcf.gz -O z -o colo829.sage.final.vcf.gz


    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA COLO829v003.sage.map.vcf.gz -O z -o COLO829v003.sage.pre.vcf.gz
    bcftools index COLO829v003.sage.pre.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA COLO829v003.sage.pre.vcf.gz -O z -o COLO829v003.sage.final.vcf.gz

    */

    private final boolean addInterimReadContexts;
    private final int minQuality;
    private final SageConfig config;
    private final GenomeRegion bounds;
    private final RefSequence refGenome;
    private final RefContextCandidates candidates;

    RefContextConsumer(boolean addInterimReadContexts, @NotNull final SageConfig config, @NotNull final GenomeRegion bounds,
            @NotNull final RefSequence refGenome, @NotNull final RefContextCandidates candidates) {
        this.bounds = bounds;
        this.refGenome = refGenome;
        this.minQuality = config.minMapQuality();
        this.candidates = candidates;
        this.addInterimReadContexts = addInterimReadContexts;

        this.config = config;
    }

    @Override
    public void accept(@NotNull final SAMRecord record) {

        if (inBounds(record)) {

            if (record.getMappingQuality() >= minQuality && !reachedDepthLimit(record)) {
                final IndexedBases refBases = refGenome.alignment(record);

                final CigarHandler handler = new CigarHandler() {
                    @Override
                    public void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
                            final int refPosition) {
                        processSnv(record, readIndex, refPosition, element.getLength(), refBases);
                    }

                    @Override
                    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
                            final int refPosition) {
                        processInsert(element, record, readIndex, refPosition, refBases);
                    }

                    @Override
                    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
                            final int refPosition) {
                        processDel(element, record, readIndex, refPosition, refBases);
                    }
                };

                CigarTraversal.traverseCigar(record, handler);
            }
        }
    }

    public void processTargeted(@NotNull final VariantHotspot hotspot, @NotNull final SAMRecord record) {

        if (inBounds(record)) {

            if (record.getMappingQuality() >= minQuality && !reachedDepthLimit(record)) {
                final IndexedBases refBases = refGenome.alignment(record);

                final CigarHandler handler = new CigarHandler() {
                    @Override
                    public void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
                            final int refPosition) {
                        processMnv(hotspot, record, readIndex, refPosition, element.getLength(), refBases);
                    }
                };

                CigarTraversal.traverseCigar(record, handler);
            }
        }
    }

    private void processInsert(@NotNull final CigarElement e, @NotNull final SAMRecord record, int readIndex, int refPosition,
            final IndexedBases refBases) {
        int refIndex = refPosition - refBases.position() + refBases.index();

        if (refPosition <= bounds.end() && refPosition >= bounds.start()) {
            final String ref = new String(refBases.bases(), refIndex, 1);
            final String alt = new String(record.getReadBases(), readIndex, e.getLength() + 1);

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (refContext != null && refContext.rawDepth() < config.maxReadDepth()) {
                int baseQuality = baseQuality(readIndex, record, alt.length());
                if (addInterimReadContexts) {
                    refContext.altRead(ref, alt, baseQuality, createInsertContext(alt, refPosition, readIndex, record, refBases));
                } else {
                    refContext.altRead(ref, alt, baseQuality);
                }
            }
        }
    }

    private void processDel(@NotNull final CigarElement e, @NotNull final SAMRecord record, int readIndex, int refPosition,
            final IndexedBases refBases) {
        int refIndex = refPosition - refBases.position() + refBases.index();

        if (refPosition <= bounds.end() && refPosition >= bounds.start()) {
            final String ref = new String(refBases.bases(), refIndex, e.getLength() + 1);
            final String alt = new String(record.getReadBases(), readIndex, 1);

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (refContext != null && refContext.rawDepth() < config.maxReadDepth()) {
                int baseQuality = baseQuality(readIndex, record, 2);
                if (addInterimReadContexts) {
                    refContext.altRead(ref, alt, baseQuality, createDelContext(ref, refPosition, readIndex, record, refBases));
                } else {
                    refContext.altRead(ref, alt, baseQuality);
                }
            }
        }
    }

    private void processSnv(@NotNull final SAMRecord record, int readBasesStartIndex, int refPositionStart, int alignmentLength,
            final IndexedBases refBases) {

        int refIndex = refPositionStart - refBases.position() + refBases.index();

        for (int i = 0; i < alignmentLength; i++) {

            int refPosition = refPositionStart + i;
            int readBaseIndex = readBasesStartIndex + i;
            int refBaseIndex = refIndex + i;

            if (!inBounds(refPosition)) {
                continue;
            }

            final byte refByte = refBases.bases()[refBaseIndex];
            final String ref = String.valueOf((char) refByte);
            final byte readByte = record.getReadBases()[readBaseIndex];

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (refContext != null && refContext.rawDepth() < config.maxReadDepth()) {
                if (readByte != refByte) {
                    int baseQuality = record.getBaseQualities()[readBaseIndex];
                    final String alt = String.valueOf((char) readByte);
                    if (addInterimReadContexts) {
                        refContext.altRead(ref, alt, baseQuality, createSNVContext(refPosition, readBaseIndex, record, refBases));
                    } else {
                        refContext.altRead(ref, alt, baseQuality);
                    }

                } else {
                    refContext.refRead();
                }
            }
        }
    }

    private void processMnv(@NotNull final VariantHotspot mnv, @NotNull final SAMRecord record, int readBasesStartIndex,
            int refPositionStart, int alignmentLength, final IndexedBases refBases) {

        final int refPositionEnd = refPositionStart + alignmentLength - 1;
        if (refPositionStart <= mnv.position() && refPositionEnd >= mnv.end()) {
            int indexOffset = (int) (mnv.position() - refPositionStart);
            int refIndex = refPositionStart - refBases.position() + refBases.index();

            int mnvRefIndex = refIndex + indexOffset;
            int mnvReadIndex = readBasesStartIndex + indexOffset;

            final String ref = new String(refBases.bases(), mnvRefIndex, mnv.ref().length());
            final String alt = new String(record.getReadBases(), mnvReadIndex, mnv.ref().length());
            if (alt.equals(mnv.alt())) {

                final RefContext refContext = candidates.refContext(record.getContig(), mnv.position());
                if (refContext != null && refContext.rawDepth() < config.maxReadDepth()) {
                    int baseQuality = baseQuality(mnvReadIndex, record, mnv.alt().length());
                    if (addInterimReadContexts) {
                        refContext.altRead(ref,
                                alt,
                                baseQuality,
                                createMNVContext((int) mnv.position(), mnvReadIndex, mnv.alt().length(), record, refBases));
                    } else {
                        refContext.altRead(ref, alt, baseQuality);
                    }
                }

            }
        }
    }

    private int baseQuality(int readIndex, SAMRecord record, int length) {
        int maxIndex = Math.min(readIndex + length, record.getBaseQualities().length) - 1;
        int quality = Integer.MAX_VALUE;
        for (int i = readIndex; i <= maxIndex; i++) {
            quality = Math.min(quality, record.getBaseQualities()[i]);
        }
        return quality;
    }

    private boolean inBounds(final SAMRecord record) {
        return record.getEnd() >= bounds.start() && record.getStart() <= bounds.end();
    }

    private boolean inBounds(final long position) {
        return position >= bounds.start() && position <= bounds.end();
    }

    private boolean reachedDepthLimit(@NotNull final SAMRecord record) {
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();

        RefContext startRefContext = candidates.refContext(bounds.chromosome(), alignmentStart);
        RefContext endRefContext = candidates.refContext(bounds.chromosome(), alignmentEnd);

        return startRefContext != null && endRefContext != null && startRefContext.rawDepth() >= config.maxReadDepth()
                && endRefContext.rawDepth() >= config.maxReadDepth();
    }

}
