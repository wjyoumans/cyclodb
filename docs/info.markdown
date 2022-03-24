---
layout: page
usemathjax: true
---

<div class="info" id="info" style="text-align:left">
</div>


<script>
const urlParams = new URLSearchParams(window.location.search);
let c = parseInt(urlParams.get('c'));

var json, id, x, text;
json = JSON.parse('{{ site.data.cyclodata.data | jsonify }}');
id = document.getElementById("info");
for (x of json) {
    if (x.conductor.plain == c) {
        text = `
            <center>
                <p>\\(K = \\mathbb{Q}(\\zeta_{${c}})\\)</p>
            <table style="width: 100%">
                <colgroup>
                    <col span="1" style="width: 30%;">
                    <col span="1" style="width: 70%;">
                </colgroup>
                <tr>
                    <td> Degree </td>
                    <td> ${x.degree.display}</td>
                </tr>
                <tr>
                    <td> Discriminant </td>
                    <td> ${x.discriminant.display} </td>
                </tr>
                <tr>
                    <td> Discriminant bits </td>
                    <td> ${x.discriminant_bits.display} </td>
                </tr>
                <tr>
                    <td> Signature </td>
                    <td> ${x.signature.display} </td>
                </tr>
                <tr>
                    <td> Galois group </td>
                    <td> ${x.galois_group.display} </td>
                </tr>
                <tr>
                    <td> Norm relation steps </td>
                    <td> ${x.norm_relation.display} </td>
                </tr>
                <tr>
                    <td> Class group </td>
                    <td> ${x.class_group.display} </td>
                </tr>
                <tr>
                    <td> Class number \\((h)\\) </td>
                    <td> ${x.h.display} </td>
                </tr>
                <tr>
                    <td> \\(h^-\\) </td>
                    <td> ${x.h_minus.display} </td>
                </tr>
                <tr>
                    <td> \\(h^+\\) </td>
                    <td> ${x.h_plus.display} </td>
                </tr>
                <tr>
                    <td> Regulator </td>
                    <td> ${x.regulator.display} </td>
                </tr>
                <tr>
                    <td> Residue </td>
                    <td> ${x.residue.display} </td>
                </tr>
            </table>
        `;
        id.innerHTML = text;
    }
}
MathJax.Hub.Queue(["Typeset",MathJax.Hub,id]);
</script>
