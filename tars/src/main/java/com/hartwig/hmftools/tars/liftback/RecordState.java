package com.hartwig.hmftools.tars.liftback;

// Whether a record reached the discriminator at all. RESOLVED means the Outcome/DecidingFeature pair
// is meaningful; the other three are record states that bypass the contest entirely and drive how the
// SAM record is rewritten (see LiftBackRecordOps).
public enum RecordState
{
    RESOLVED,
    UNMAPPED,
    LIFT_FAILED,
    SUPPLEMENTARY
}
