---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
usemathjax: true
---

<script src="https://code.jquery.com/jquery-3.2.1.min.js"></script>
<script src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>

<button type="button" class="advanced-search-button" id="advanced-search-button">
    Advanced search
</button>
<div class="content" id="advanced-search">
    <p>
        Advanced search goes here. Ideas: galois group/class group iscyclic checkbox, galois 
        group/class group prime factors, upper/lower bounds where applicable
    </p>
</div>

<script>
let coll = document.getElementById("advanced-search-button");
let content = document.getElementById("advanced-search");
coll.addEventListener("click", function() {    
    this.classList.toggle("active");
    if (content.style.display === "block") {
      content.style.display = "none";
    } else {
      content.style.display = "block";
    }
});
</script>


<table id="table" class="display nowrap" width="100%">
<thead>
  <tr>
    <th>Conductor</th>
    <th>Degree</th>
    <th>Norm relation</th>
    <th>Signature</th>
    <th>Galois group</th>
    <th>Discriminant</th>
    <th>Discriminant bits</th>
    <th>Class number (\(h\))</th>
    <th>\(h^-\)</th>
    <th>\(h^+\)</th>
    <th>Class group</th>
    <th>Regulator</th>
    <th>Residue</th>
    <th>Polynomial</th>
    <th>Precision</th>
  </tr>
</thead>
<tfoot>
  <tr>
    <th>Conductor</th>
    <th>Degree</th>
    <th>Norm relation</th>
    <th>Signature</th>
    <th>Galois group</th>
    <th>Discriminant</th>
    <th>Discriminant bits</th>
    <th>Class number (\(h\))</th>
    <th>\(h^-\)</th>
    <th>\(h^+\)</th>
    <th>Class group</th>
    <th>Regulator</th>
    <th>Residue</th>
    <th>Polynomial</th>
    <th>Precision</th>
  </tr>
</tfoot>
</table>

<script>
$(document).ready(function() {
    var json = JSON.parse('{{ site.data.cyclodata.data | jsonify }}');
    var id = document.getElementById("table");
    var table = $('#table').DataTable({
        data: json,
        pageLength: 10,
        columnDefs: [{
            targets: [1],
            orderData: [1, 0],
        },
        {
            targets: [2],
            orderData: [2, 1, 0],
        },
        {
            targets: [3, 4],
            orderData: [1, 0],
        },
        {
            targets: [7, 10],
            orderData: [7, 1, 0],
        },
        {
            targets: [8],
            orderData: [8, 1, 0],
        },
        {
            targets: [9],
            orderData: [9, 1, 0],
        }
        ],
        columns: [
            { 
                data: "conductor", 
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "degree", 
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "norm_relation", 
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "signature", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "galois_group", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "discriminant", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "discriminant_bits", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "h", 
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "h_minus",
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "h_plus", 
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "class_group", 
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "regulator",
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "residue",
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "polynomial", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
            { 
                data: "precision",
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"}
            },
        ],
        drawCallback: function(settings) {
            MathJax.Hub.Queue(["Typeset",MathJax.Hub,id]);
        },
        headerCallback: function(settings) {
            MathJax.Hub.Queue(["Typeset",MathJax.Hub,id]);
        }
    });
    $('#table tbody').on('click', 'tr', function () {
        let data = table.row(this).data();
        window.location = `/info?c=${data.conductor.plain}`
    });
});

</script>
