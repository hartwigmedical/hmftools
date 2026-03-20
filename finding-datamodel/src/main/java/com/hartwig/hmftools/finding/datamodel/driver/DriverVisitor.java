package com.hartwig.hmftools.finding.datamodel.driver;

import com.hartwig.hmftools.finding.datamodel.ChromosomeArmCopyNumber;
import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.Virus;

/**
 * Visitor interface for processing the concrete subtypes of {@link Driver}.
 * <p>
 * Because {@code Driver} is a sealed type hierarchy (SmallVariant, GainDeletion, Fusion,
 * Disruption, Virus, ChromosomeArmCopyNumber), this visitor lets callers handle each subtype
 * without casting. Each {@code Driver} implementation calls back the matching {@code visit}
 * overload via its {@code accept} method.
 * <p>
 * Usage: implement this interface and pass the implementation to {@link Driver#accept}.
 * <pre>{@code
 * class MyVisitor implements DriverVisitor {
 *     public void visit(SmallVariant v)              { ... }
 *     public void visit(GainDeletion g)              { ... }
 *     public void visit(Fusion f)                    { ... }
 *     public void visit(Disruption d)                { ... }
 *     public void visit(Virus v)                     { ... }
 *     public void visit(ChromosomeArmCopyNumber c)   { ... }
 * }
 * driver.accept(new MyVisitor());
 * }</pre>
 *
 * @see DriverGenesGetter for a concrete example
 */
public interface DriverVisitor
{
    void visit(SmallVariant smallVariant);
    void visit(GainDeletion gainDeletion);
    void visit(Fusion fusion);
    void visit(Disruption disruption);
    void visit(Virus virus);
    void visit(ChromosomeArmCopyNumber chromosomeArmCopyNumber);
}
