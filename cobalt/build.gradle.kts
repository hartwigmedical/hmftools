
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - COBALT"

dependencies {
    implementation(project(":hmf-common"))
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.cobalt.CobaltApplication")
}
