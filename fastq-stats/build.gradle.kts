
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - FASTQ Stats"
version = "1.2"

dependencies {
    implementation(libs.annotations)
    implementation(libs.log4j.core)
    implementation(libs.google.guava)
    implementation(libs.commons.cli)
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.fastqstats.FastqStatsRunner")
}