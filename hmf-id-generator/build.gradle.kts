
plugins {
    kotlin("jvm")
    id("com.hartwig.kotlin-conventions")
}

description = "HMF Tools - HMF ID Generator"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    implementation(libs.bcprov.jdk15on)
    implementation(libs.kotlin.reflect)
    implementation(libs.commons.csv)
    implementation(libs.commons.cli)
    implementation(libs.selenium.java)
    testImplementation(libs.junit)
    testImplementation(libs.kotest.runner.junit5.jvm)
    testImplementation(libs.kotest.assertions.core.jvm)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.idgenerator.HmfIdApplicationKt")
}
