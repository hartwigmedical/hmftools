package com.hartwig.hmftools.civic.data;

public abstract class CivicVariantKey {
    public abstract String name();

    public abstract String id();

    public abstract int acceptedEvidenceItems();

    public abstract int rejectedEvidenceItems();

    public abstract int submittedEvidenceItems();

}
