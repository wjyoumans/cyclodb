---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
usemathjax: true
---


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
        searchBuilder: {
            greyscale: true,
            columns: [0,1,2,5,6,7,8,9,11,12],
            conditions:{
                num:{
                    'MultipleOf': {
                        conditionName: 'Multiple Of',
                        init: function (that, fn, preDefined = null) {
                            // Declare the input element and set the listener to trigger searching
                            var el =  $('<input/>').on('input', function() { fn(that, this) });
     
                            // Add mechanism to apply preDefined values that may be passed in
                            if (preDefined !== null) {
                                $(el).val(preDefined[0]);
                            }
     
                            return el;
                        },
                        inputValue: function (el) {
                            // Return the value within the input element
                            return $(el[0]).val();
                        },
                        isInputValid: function (el, that) {
                            // If there is text in the input element then it is valid for searching
                            return $(el[0]).val().length !== 0;
                        },
                        search: function (value, comparison) {
                            // Use the modulo (%) operator to check that there is no remainder
                            return value%comparison === 0;
                        }
                    }
                }
            }
        },
        dom: 'Qlfrtip',
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
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "degree", 
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "norm_relation", 
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "signature", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "galois_group", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "discriminant",
                render: {"filter": "filter", "display": "display", "_": "plain"},
                type: "num",
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "discriminant_bits", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "h", 
                render: {"filter": "filter", "display": "display", "_": "plain"},
                type: "num",
                searchBuilderTitle: "h",
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "h_minus",
                visible: false,
                render: {"filter": "filter", "display": "display", "_": "plain"},
                type: "num",
                searchBuilderTitle: "h-",
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "h_plus", 
                render: {"filter": "filter", "display": "display", "_": "plain"},
                type: "num",
                searchBuilderTitle: "h+",
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "class_group", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "regulator",
                render: {"filter": "filter", "display": "display", "_": "plain"},
                type: "num",
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "residue",
                visible: false,
                render: {"filter": "filter", "display": "display", "_": "plain"},
                type: "num",
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "polynomial", 
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
            },
            { 
                data: "precision",
                visible: false,
                searchable: false,
                render: {"filter": "filter", "display": "display", "_": "plain"},
                searchBuilder: {
                    orthogonal: {
                        display: "filter",
                    }
                },
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
        window.location = '{{ "/info" | relative_url }}' + `?c=${data.conductor.plain}`
    });
    $("#deg_min").keyup(function() {
        $.fn.dataTable.ext.search.push(
            function(settings, data, dataIndex) {
                return data[1] >= parseInt($("#deg_min").val());
            }
        );
        table.draw();
    });
    $("#deg_max").keyup(function() {
        $.fn.dataTable.ext.search.push(
            function(settings, data, dataIndex) {
                return data[1] <= parseInt($("#deg_max").val());
            }
        );
        table.draw();
    });
    $("#disc_min").keyup(function() {
        $.fn.dataTable.ext.search.push(
            function(settings, data, dataIndex) {
                return data[6] >= parseInt($("#disc_min").val());
            }
        );
        table.draw();
    });
    $("#disc_max").keyup(function() {
        $.fn.dataTable.ext.search.push(
            function(settings, data, dataIndex) {
                return data[6] <= parseInt($("#disc_max").val());
            }
        );
        table.draw();
    });
    $("#h_min").keyup(function() {
        $.fn.dataTable.ext.search.push(
            function(settings, data, dataIndex) {
                return data[7] >= parseInt($("#h_min").val());
            }
        );
        table.draw();
    });
    $("#h_max").keyup(function() {
        $.fn.dataTable.ext.search.push(
            function(settings, data, dataIndex) {
                return data[7] <= parseInt($("#h_max").val());
            }
        );
        table.draw();
    });
    $("#reg_min").keyup(function() {
        $.fn.dataTable.ext.search.push(
            function(settings, data, dataIndex) {
                return data[11] >= parseInt($("#reg_min").val());
            }
        );
        table.draw();
    });
    $("#reg_max").keyup(function() {
        $.fn.dataTable.ext.search.push(
            function(settings, data, dataIndex) {
                return data[11] <= parseInt($("#reg_max").val());
            }
        );
        table.draw();
    });
});
</script>
