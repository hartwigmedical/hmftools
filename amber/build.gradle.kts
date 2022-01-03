
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - AMBER"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.commons.cli)
    implementation(libs.jcommander)
    implementation(libs.tablesaw.core)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.amber.AmberApplication")
}
