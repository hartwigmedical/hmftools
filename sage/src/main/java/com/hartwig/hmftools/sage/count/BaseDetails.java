package com.hartwig.hmftools.sage.count;

import java.util.List;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.ModifiableVariantHotspotEvidence;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidence;

import org.jetbrains.annotations.NotNull;

public class BaseDetails implements Comparable<BaseDetails> {

    private final String contig;
    private final long position;
    private int refSupport;
    private int refQuality;
    private final List<ModifiableVariantHotspotEvidence> evidenceList = Lists.newArrayList();

    public BaseDetails(final String contig, final long position) {
        this.contig = contig;
        this.position = position;
    }

    public long position() {
        return position;
    }

    public boolean isEmpty() {
        return evidenceList.isEmpty();
    }

    public List<VariantHotspotEvidence> evidence() {
        List<VariantHotspotEvidence> result = Lists.newArrayList();
        for (ModifiableVariantHotspotEvidence evidence : evidenceList) {
            result.add(update(evidence));
        }

        return result;
    }

    private ModifiableVariantHotspotEvidence update(ModifiableVariantHotspotEvidence evidence) {
        return evidence.setRefQuality(refQuality).setRefSupport(refSupport).setReadDepth(refSupport + evidence.altSupport());
    }

    @NotNull
    public ModifiableVariantHotspotEvidence selectOrCreate(@NotNull final String ref, @NotNull final String alt) {
        for (ModifiableVariantHotspotEvidence evidence : evidenceList) {
            if (evidence.ref().equals(ref) && evidence.alt().equals(alt)) {
                return update(evidence);
            }
        }

        ModifiableVariantHotspotEvidence newEvidence = ModifiableVariantHotspotEvidence.create()
                .setChromosome(contig)
                .setPosition(position)
                .setRef(ref)
                .setAlt(alt)
                .setAltSupport(0)
                .setAltQuality(0)
                .setIndelSupport(0);

        evidenceList.add(update(newEvidence));
        return newEvidence;
    }

    public void incrementRefSupport() {
        this.refSupport++;
    }

    public void incrementRefQuality(final int quality) {
        this.refQuality += quality;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final BaseDetails that = (BaseDetails) o;
        return position == that.position;
    }

    @Override
    public int hashCode() {
        return Objects.hash(position);
    }

    @Override
    public int compareTo(@NotNull final BaseDetails o) {
        return Long.compare(position, o.position);
    }
}
