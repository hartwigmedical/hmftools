
plugins {
    id("com.hartwig.java-conventions")
    id("com.hartwig.tests-jar") // create a jar for tests
}

description = "HMF Tools - Common"

dependencies {
    api(libs.annotations)
    api(libs.commons.cli)
    api(libs.log4j.core)
    api(libs.google.guava)
    api(libs.google.gson)
    api(libs.htsjdk)
    api(libs.commons.math3)
    api(libs.commons.lang3)
    api(libs.jcommander)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    compileOnly(libs.immutables.gson)
    testImplementation(libs.junit)
}
