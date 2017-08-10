package com.hartwig.hmftools.civic.data;

import java.util.List;

public abstract class CivicVariant {
    public abstract String id();

    public abstract String gene();

    public abstract String name();

    public abstract String description();

    public abstract String chromosome();

    public abstract String start();

    public abstract String stop();

    public abstract String representativeTranscript();

    public abstract String referenceBuild();

    public abstract List<CivicEvidenceItem> evidenceItems();

}
