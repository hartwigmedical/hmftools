package com.hartwig.hmftools.finding.datamodel;

import java.lang.annotation.ElementType;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@io.soabase.recordbuilder.core.RecordBuilder.Template(options = @io.soabase.recordbuilder.core.RecordBuilder.Options(
        nullableAnnotationClass = "org.jspecify.annotations.Nullable",
        defaultNotNull = true,
        useImmutableCollections = true
        //interpretNotNulls = true
))
@Retention(RetentionPolicy.SOURCE)
@Target(ElementType.TYPE)
@Inherited
public @interface RecordBuilder {
}